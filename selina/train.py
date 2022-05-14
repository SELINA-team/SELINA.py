import torch
from torch.autograd import Function
import torch.nn as nn
import glob
import ntpath
import datatable as dt
import pandas as pd
import numpy as np
from functools import reduce
from torch.utils.data import Dataset
from tqdm import tqdm
from torch.utils.data import DataLoader
import torch.nn.functional as F
import argparse as ap
import pickle
from imblearn.over_sampling import SMOTE
from collections import Counter
import os

params_train = [0.0001, 50, 128]
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')


def train_parser(subparsers):
    train = subparsers.add_parser(
        "train", help="Pre-training using provided data based on MADA.")
    train.add_argument('--path-in',
                       dest='path_in',
                       type=str,
                       required=True,
                       help='File path of training datasets.')
    train.add_argument('--path-out',
                       dest='path_out',
                       type=str,
                       required=True,
                       help='File path of the output model.')
    train.add_argument('--outprefix',
                       type=str,
                       required=False,
                       help='Prefix of the output files. DEFAULT: pre-trained',
                       default='pre-trained')


def read_expr(path):
    expr = dt.fread(path, header=True, sep='\t', nthreads=6)
    expr = expr.to_pandas()
    expr.index = expr.loc[:, 'Gene']
    del expr['Gene']
    expr = expr.astype(float)
    return expr


def label2dic(label):
    label_set = list(set(label))
    dic = {}
    for i in range(len(label_set)):
        dic[label_set[i]] = i
    return dic


def preprocessing(path_in):
    samples = sorted(glob.glob(path_in + '/*_expr.txt'))
    metas = sorted(glob.glob(path_in + '/*_meta.txt'))
    train_sets = []
    celltypes = []
    platforms = []
    genes = []
    print('Loading data')
    for i in range(len(samples)):
        train_sets.append(read_expr(samples[i]))
        meta = pd.read_csv(metas[i], sep='\t', header=0)
        celltype = meta['Celltype'].to_list()
        platform = meta['Platform'].to_list()
        celltypes.append(celltype)
        platforms.append(platform)
        genes.append(train_sets[i].index.to_list())
    genes = reduce(np.intersect1d, genes)
    ct_freqs = Counter([i for item in celltypes for i in item])
    max_n = max(ct_freqs.values())
    rct_freqs = {}
    if max_n < 500:
        sample_n = 100
    elif max_n < 1000:
        sample_n = 500
    else:
        sample_n = 1000
    for ct, freq in ct_freqs.items():
        if freq <= sample_n:
            rct_freqs[ct] = freq

    for i in range(len(samples)):
        sample_ct_freq = {}
        ct_freq = Counter(celltypes[i])
        if len(ct_freq) > 1:
            for ct, freq in rct_freqs.items():
                if (ct in ct_freq.keys()) & (ct_freq[ct] >= 6):
                    sample_ct_freq[ct] = round(sample_n * ct_freq[ct] / freq)
            smo = SMOTE(sampling_strategy=sample_ct_freq, random_state=1)
            train_sets[i], celltypes[i] = smo.fit_resample(
                train_sets[i].T, celltypes[i])
            train_sets[i] = train_sets[i].T
            platforms[i] = np.unique(
                platforms[i]).tolist() * train_sets[i].shape[1]
    platforms = [i for item in platforms for i in item]
    celltypes = [i for item in celltypes for i in item]
    for i in range(len(samples)):
        train_sets[i] = train_sets[i].loc[genes, ]
        train_sets[i] = np.divide(train_sets[i], np.sum(train_sets[i],
                                                        axis=0)) * 10000
        train_sets[i] = np.log2(train_sets[i] + 1)
    train_data = pd.concat(train_sets, axis=1)
    return train_data, celltypes, platforms, genes


# def preprocessing(path_in):
#     samples = sorted(glob.glob(path_in + '/*_expr.txt'))
#     metas = sorted(glob.glob(path_in + '/*_meta.txt'))
#     train_sets = []
#     celltypes = []
#     platforms = []
#     genes = []
#     print('Loading data')
#     for i in range(len(samples)):
#         train_sets.append(read_expr(samples[i]))
#         genes.append(train_sets[i].index.to_list())
#         meta = pd.read_csv(metas[i], sep='\t', header=0)
#         celltypes = celltypes + meta['Celltype'].to_list()
#         platforms = platforms + meta['Platform'].to_list()
#     genes = reduce(np.intersect1d, genes)
#     for i in range(len(samples)):
#         train_sets[i] = train_sets[i].loc[genes, ]
#         train_sets[i] = np.divide(train_sets[i], np.sum(train_sets[i],
#                                                         axis=0)) * 10000
#         train_sets[i] = np.log2(train_sets[i] + 1)
#     train_data = pd.concat(train_sets, axis=1)
#     return train_data, celltypes, platforms, genes


class Datasets(Dataset):
    def __init__(self, data, celltypes, platforms, ct_dic, plat_dic):
        class_labels = [ct_dic[i] for i in celltypes]
        domain_labels = [plat_dic[i] for i in platforms]
        self.class_labels = torch.as_tensor(class_labels)
        self.domain_labels = torch.as_tensor(domain_labels)
        self.expr = data.values

    def __getitem__(self, index):
        return torch.as_tensor(
            self.expr[:, index]
        ), self.class_labels[index], self.domain_labels[index]

    def __len__(self):
        return len(self.class_labels)


class GRL(Function):
    @staticmethod
    def forward(ctx, x, alpha):
        ctx.alpha = alpha
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        output = grad_output.neg() * ctx.alpha
        return output, None


class MADA(nn.Module):
    def __init__(self, nfeatures, nct, nplat):
        super(MADA, self).__init__()
        self.nct = nct
        self.feature = nn.Sequential(
            nn.Linear(in_features=nfeatures, out_features=100), nn.ReLU(),
            nn.Dropout())
        self.class_classifier = nn.Sequential(
            nn.Linear(in_features=100, out_features=50), nn.ReLU(),
            nn.Dropout(), nn.Linear(in_features=50, out_features=nct))
        self.domain_classifier = nn.ModuleList([
            nn.Sequential(nn.Linear(in_features=100, out_features=25),
                          nn.ReLU(),
                          nn.Linear(in_features=25, out_features=nplat))
            for _ in range(nct)
        ])

    def forward(self, input_data, alpha, nct):
        features = self.feature(input_data)
        class_logits = self.class_classifier(features)
        class_predictions = F.softmax(class_logits, dim=1)
        reverse_features = GRL.apply(features, alpha)
        domain_logits = []
        for class_idx in range(nct):
            wrf = class_predictions[:,
                                    class_idx].unsqueeze(1) * reverse_features
            domain_logits.append(self.domain_classifier[class_idx](wrf))
        return class_logits, domain_logits


def train(train_data, params, celltypes, platforms, nfeatures, nct, nplat,
          ct_dic, plat_dic, device):
    network = MADA(nfeatures, nct, nplat).train()
    lr = params[0]
    n_epoch = params[1]
    batch_size = params[2]
    train_data = Datasets(train_data, celltypes, platforms, ct_dic, plat_dic)
    optimizer = torch.optim.Adam(network.parameters(), lr=lr)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    network = network.to(device)
    loss_class = loss_class.to(device)
    loss_domain = loss_domain.to(device)
    train_loader = DataLoader(dataset=train_data,
                              batch_size=batch_size,
                              shuffle=True,
                              drop_last=True)

    len_train_loader = len(train_loader)
    print('Begin training')
    for epoch in tqdm(range(n_epoch)):
        loader_iter = iter(train_loader)
        for i in range(len_train_loader):
            p = float(i +
                      epoch * len_train_loader) / n_epoch / len_train_loader
            alpha = 2. / (1. + np.exp(-10 * p)) - 1
            expr, class_label, domain_label = loader_iter.next()
            expr = expr.to(device)
            expr = expr.float()
            class_label = class_label.to(device)
            domain_label = domain_label.to(device)
            class_output, domain_output = network(input_data=expr,
                                                  alpha=alpha,
                                                  nct=nct)
            err_class = loss_class(class_output, class_label)
            err_domain = [
                loss_domain(domain_output[class_idx], domain_label)
                for class_idx in range(nct)
            ]
            loss_total = (1 -
                          alpha) * sum(err_domain) / nct + alpha * err_class
            optimizer.zero_grad()
            loss_total.backward()
            optimizer.step()
    print('Finish Training')
    return network


def train_model(path_in, path_out, outprefix):
    train_data, celltypes, platforms, genes = preprocessing(path_in)
    ct_dic, plat_dic = label2dic(celltypes), label2dic(platforms)
    nfeatures, nct, nplat = train_data.shape[0], len(ct_dic), len(plat_dic)
    network = train(train_data, params_train, celltypes, platforms, nfeatures,
                    nct, nplat, ct_dic, plat_dic, device)

    try:
        os.makedirs(path_out)
    except OSError:
        pass

    torch.save(network, path_out + '/' + outprefix + '_params.pt')
    model_meta = {'genes': genes, 'celltypes': ct_dic}
    with open(path_out + '/' + outprefix + '_meta.pkl', 'wb') as f:
        pickle.dump(model_meta, f)
    print('All done')