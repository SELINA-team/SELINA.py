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
    train.add_argument(
        '--disease',
        action="store_true",
        help=
        'This flag should be used when the data is in some disease condition')

params_train = [0.0001, 50, 128]
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

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

def delete_multiple_element(list_object, indices):
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)
    return list_object

def preprocessing(path_in, disease):
    samples = sorted(glob.glob(path_in + '/*_expr.txt'))
    metas = sorted(glob.glob(path_in + '/*_meta.txt'))
    train_sets = []
    celltypes = []
    if disease:
        diseases = []
    platforms = []
    genes = []
    print('Loading data')
    
    for i in range(len(samples)):
        train_sets.append(read_expr(samples[i]))
        meta = pd.read_csv(metas[i], sep='\t', header=0)
        celltypes.append(meta['Celltype'].to_numpy())
        platforms.append(meta['Platform'].to_numpy())
        genes.append(train_sets[i].index.to_list())
        if disease:
            diseases.append(meta['Disease'].to_numpy())


    genes = reduce(np.intersect1d, genes)

    length = len(train_sets)
    drop_index = []
    for i in range(length):
        sample_plat = np.unique(platforms[i])
        if len(sample_plat)>=2 :
            drop_index.append(i)
            for plat in sample_plat:
                index = np.where(platforms[i]==plat)[0]
                train_sets.append(train_sets[i].iloc[:,index])
                celltypes.append(celltypes[i][index])
                if disease:
                    diseases.append(diseases[i][index])
                platforms.append(platforms[i][index])
    train_sets = delete_multiple_element(train_sets,drop_index)
    celltypes = delete_multiple_element(celltypes,drop_index)
    if disease:
        diseases = delete_multiple_element(diseases,drop_index)
    platforms =delete_multiple_element(platforms,drop_index)
    
    if disease:
        length = len(train_sets)
        drop_index = []
        for i in range(length):
            disease_type = np.unique(diseases[i])
            if len(disease_type) >=2 :
                drop_index.append(i)
                for dis in disease_type:
                    index = np.where(np.array(diseases[i])==dis)[0]
                    train_sets.append(train_sets[i].iloc[:,index])
                    celltypes.append(celltypes[i][index])
                    diseases.append(diseases[i][index])
                    platforms.append(platforms[i][index])
        train_sets = delete_multiple_element(train_sets,drop_index)
        celltypes = delete_multiple_element(celltypes,drop_index)
        diseases = delete_multiple_element(diseases,drop_index)
        platforms =delete_multiple_element(platforms,drop_index)

    ct_freqs = Counter([i for item in celltypes for i in item])
    max_n = max(ct_freqs.values())
    rct_freqs = {}
    if  (max_n < 500) & (max_n) >100:
        sample_n = 100
    elif max_n < 1000:
        sample_n = 500
    else:
        sample_n = 1000
    for ct,freq in ct_freqs.items():
        if freq <= sample_n:
            rct_freqs[ct] = freq

    for i in range(len(train_sets)):                                            
        sample_ct_freq = {}
        ct_freq = Counter(celltypes[i])
        if len(ct_freq)>1:
            for ct,freq in rct_freqs.items():
                if (ct in ct_freq.keys()) & (ct_freq[ct] >= 4):
                    sample_ct_freq[ct] = round(sample_n * ct_freq[ct]/freq)
            smo = SMOTE(sampling_strategy = sample_ct_freq,random_state=1, k_neighbors=3)
            train_sets[i],celltypes[i] = smo.fit_resample(train_sets[i].T,celltypes[i])
            train_sets[i] = train_sets[i].T
            platforms[i] = np.unique(platforms[i]).tolist() * train_sets[i].shape[1]
            if disease:
                diseases[i] = np.unique(diseases[i]).tolist() * train_sets[i].shape[1]
    
    celltypes = [i for item in celltypes for i in item]
    if disease:
        diseases = [i for item in diseases for i in item]
    platforms = [i for item in platforms for i in item]
    
    for i in range(len(train_sets)):
        train_sets[i] = train_sets[i].loc[genes, ]
        train_sets[i] = np.divide(train_sets[i], np.sum(train_sets[i],
                                                        axis=0)) * 10000
        train_sets[i] = np.log2(train_sets[i] + 1)
    train_data = pd.concat(train_sets, axis=1)
    if disease:
        return train_data, celltypes, diseases, platforms, genes
    else:
        return train_data, celltypes, platforms, genes

class Datasets(Dataset):
    def __init__(self, data, celltypes, diseases, platforms, ct_dic, disease_dic, plat_dic, disease):
        class_labels = [ct_dic[i] for i in celltypes]
        domain_labels = [plat_dic[i] for i in platforms]
        self.class_labels = torch.as_tensor(class_labels)
        self.domain_labels = torch.as_tensor(domain_labels)
        if disease:
            self.disease_flag = True
            disease_labels =  [disease_dic[i] for i in diseases]
            self.disease_labels = torch.as_tensor(disease_labels)
        else:
            self.disease_flag = False
        self.expr = data.values

    def __getitem__(self, index):
        if self.disease_flag:
            return torch.as_tensor(
            self.expr[:, index]
        ), self.class_labels[index], self.disease_labels[index],  self.domain_labels[index]
        else:
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
    def __init__(self, nfeatures, nct, ndis, nplat, disease):
        super(MADA, self).__init__()
        if disease:
            self.disease_flag = True
            self.nct = nct
            self.feature = nn.Sequential(
                nn.Linear(in_features=nfeatures, out_features=500), nn.ReLU(),
                nn.Dropout())
            self.disease_classifier = nn.Sequential(
                nn.Linear(in_features=500, out_features=50), nn.ReLU(),
                nn.Linear(in_features=50, out_features=ndis))
            self.class_classifier = nn.Sequential(
                nn.Linear(in_features=500, out_features=50), nn.ReLU(),
                nn.Dropout(), nn.Linear(in_features=50, out_features=nct))
            self.domain_classifier = nn.ModuleList([
                nn.Sequential(nn.Linear(in_features=500, out_features=25),
                            nn.ReLU(),
                            nn.Linear(in_features=25, out_features=nplat))
                for _ in range(nct+ndis)
            ])
        else:
            self.disease_flag = False
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
    def forward(self, input_data, alpha, nct, ndis):
        features = self.feature(input_data)
        class_logits = self.class_classifier(features)
        class_predictions = F.softmax(class_logits, dim=1)
        reverse_features = GRL.apply(features, alpha)
        domain_logits = []
        for class_idx in range(nct):
            wrf = class_predictions[:,
                                    class_idx].unsqueeze(1) * reverse_features
            domain_logits.append(self.domain_classifier[class_idx](wrf))
        if self.disease_flag:
            disease_logits =  self.disease_classifier(features)
            disease_preditcions =  F.softmax(disease_logits, dim=1)
            for dis_idx in range(ndis):
                wrf = disease_preditcions[:,
                                        dis_idx].unsqueeze(1) * reverse_features
                domain_logits.append(self.domain_classifier[dis_idx](wrf))
            return class_logits, disease_logits, domain_logits
        else:
            return class_logits, domain_logits


def train(train_data, params, celltypes, diseases, platforms, nfeatures, nct, ndis,nplat, ct_dic, disease_dic, plat_dic, device, disease):
    if disease:
        network = MADA(nfeatures, nct, ndis, nplat, disease).train()
        train_data = Datasets(train_data, celltypes, diseases, platforms, ct_dic, disease_dic, plat_dic, disease)
    else:
        network = MADA(nfeatures, nct, None, nplat, disease).train()
        train_data = Datasets(train_data, celltypes, None, platforms, ct_dic, None, plat_dic, disease)
        loss_disease = nn.CrossEntropyLoss()
        loss_disease = loss_disease.to(device)
    lr = params[0]
    n_epoch = params[1]
    batch_size = params[2]
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
            alpha = 2. / (1. + np.exp(-10* p)) - 1
            if disease:
                expr, class_label, disease_label, domain_label = loader_iter.next()
                expr = expr.to(device)
                expr = expr.float()
                class_label = class_label.to(device)
                disease_label = disease_label.to(device)
                domain_label = domain_label.to(device)
                class_output, disease_output, domain_output = network(input_data=expr,
                                                    alpha=alpha,
                                                    nct=nct,
                                                    ndis=ndis)
                err_class = loss_class(class_output, class_label)
                err_disease =  loss_class(disease_output, disease_label)
                err_domain_ct = [
                    loss_domain(domain_output[class_idx], domain_label)
                    for class_idx in range(nct)
                ]
                err_domain_dis = [
                    loss_domain(domain_output[dis_idx], domain_label)
                    for dis_idx in range(nct,(nct+ndis))
                ]
                err_domain = err_domain_ct+err_domain_dis
                loss_total = (1 - alpha)/2 * (1- alpha) * (sum(err_domain_dis) / ndis) +  (alpha * 2) * (1 - alpha) * (sum(err_domain_ct) / nct) + alpha * (err_class +  err_disease)
            else:
                expr, class_label, domain_label = loader_iter.next()
                expr = expr.to(device)
                expr = expr.float()
                class_label = class_label.to(device)
                domain_label = domain_label.to(device)
                class_output, domain_output = network(input_data=expr,
                                                    alpha=alpha,
                                                    nct=nct,
                                                    ndis=None)
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


def train_model(path_in, path_out, disease, outprefix):
    if disease:
        train_data, celltypes, diseases, platforms, genes = preprocessing(path_in, disease)
        ct_dic, disease_dic, plat_dic = label2dic(celltypes), label2dic(diseases), label2dic(platforms)
        nfeatures, nct, ndis, nplat = train_data.shape[0], len(ct_dic),len(disease_dic), len(plat_dic)
        network = train(train_data, params_train, celltypes, diseases, platforms, nfeatures, nct, ndis, nplat, ct_dic, disease_dic, plat_dic, device, disease)

        try:
            os.makedirs(path_out)
        except OSError:
            pass

        torch.save(network, path_out + '/' + outprefix + '_params.pt')
        model_meta = {'genes': genes, 'cellsources': disease_dic, 'celltypes': ct_dic}
        with open(path_out + '/' + outprefix + '_meta.pkl', 'wb') as f:
            pickle.dump(model_meta, f)
        print('All done')
    else:
        train_data, celltypes, platforms, genes = preprocessing(path_in, disease)
        ct_dic, plat_dic = label2dic(celltypes), label2dic(platforms)
        nfeatures, nct, nplat = train_data.shape[0], len(ct_dic), len(plat_dic)
        network = train(train_data, params_train, celltypes, None, platforms, nfeatures, nct, None, nplat, ct_dic, None, plat_dic, device, disease)

        try:
            os.makedirs(path_out)
        except OSError:
            pass

        torch.save(network, path_out + '/' + outprefix + '_params.pt')
        model_meta = {'genes': genes, 'celltypes': ct_dic}
        with open(path_out + '/' + outprefix + '_meta.pkl', 'wb') as f:
            pickle.dump(model_meta, f)
        print('All done')