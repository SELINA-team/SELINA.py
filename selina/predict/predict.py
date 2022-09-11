from email.policy import default
import torch
from torch.autograd import Function
import torch.nn as nn
import datatable as dt
import pandas as pd
import numpy as np
from functools import reduce
from torch.utils.data import Dataset
from tqdm import tqdm
from torch.utils.data import DataLoader
import torch.nn.functional as F
import argparse as ap
from selina.train import read_expr, GRL, MADA
import pickle
from functools import reduce
import os

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
params_tune1 = [0.0005, 50, 128]
params_tune2 = [0.0001, 10, 128]


def predict_parser(subparsers):
    workflow = subparsers.add_parser(
        'predict', help='Fine-tuning and predict for the query data. ')
    group_input = workflow.add_argument_group('Arguments for input')
    group_input.add_argument('--query-expr',
                             dest='query_expr',
                             type=str,
                             required=True,
                             help='File path of the query data matrix.')
    group_input.add_argument('--model',
                             type=str,
                             required=True,
                             help='File path of the pre-trained model.')
    group_input.add_argument('--seurat',
                             type=str,
                             required=True,
                             help='File path of the seurat object.')
    group_input.add_argument(
        '--disease',
        action="store_true",
        help=
        'This flag should be used when the data is in some disease condition')
    group_cutoff = workflow.add_argument_group(
        'Cutoff for downstream analysis')
    group_cutoff.add_argument(
        '--prob-cutoff',
        dest='prob_cutoff',
        type=float,
        required=False,
        default=0.9,
        help='Cutoff for prediction probability. DEFAULT: 0.9.')

    group_output = workflow.add_argument_group('Output Arguments')
    group_output.add_argument('--path-out',
                              dest='path_out',
                              type=str,
                              required=True,
                              help='File path of the output files.')
    group_output.add_argument(
        '--outprefix',
        type=str,
        required=False,
        default='query',
        help='Prefix of the output files. DEFAULT: query')


def merge(genes, query_expr):
    model_expr = pd.DataFrame(np.random.randn(len(genes), 1))
    model_expr.columns = ['genes']
    model_expr.index = genes
    query_expr = pd.DataFrame.join(model_expr, query_expr).drop('genes',
                                                                axis=1)
    query_expr = query_expr.fillna(0)
    query_expr = np.divide(query_expr, np.sum(query_expr, axis=0)) * 10000
    query_expr = np.log2(query_expr + 1)
    return query_expr


class Datasets(Dataset):
    def __init__(self, data):
        self.expr = data.values

    def __getitem__(self, index):
        return torch.as_tensor(self.expr[:, index])

    def __len__(self):
        return self.expr.shape[1]


class Autoencoder(nn.Module):
    def __init__(self, network, nfeature, nct):
        super(Autoencoder, self).__init__()
        encoder = list(network.feature.children()) + list(
            network.class_classifier.children())
        encoder_index = [0, 1, 3, 4, 6]
        self.encoder = nn.Sequential(*[encoder[i] for i in encoder_index],
                                     nn.ReLU())
        self.decoder = nn.Sequential(
            nn.Linear(in_features=nct, out_features=50), nn.ReLU(),
            nn.Linear(in_features=50, out_features=100), nn.ReLU(),
            nn.Linear(in_features=100, out_features=nfeature))

    def forward(self, input_data):
        output = self.decoder(self.encoder(input_data))
        return (output)


class Normal_Classifier(nn.Module):
    def __init__(self, network):
        super(Normal_Classifier, self).__init__()
        self.classifier = nn.Sequential(*list(network.encoder.children())[:-1],
                                        nn.Softmax(dim=1))

    def forward(self, input_data):
        output = self.classifier(input_data)
        return (output)

class Disease_Classifier(nn.Module):
    def __init__(self,network):
        super(Disease_Classifier,self).__init__()
        self.feature = nn.Sequential(*[list(network.feature.children())[i] for i in  [0,1]])
        self.celltype =  nn.Sequential(*[list(network.class_classifier.children())[i] for i in  [0,1,3]], nn.Softmax(dim=1))
        self.disease =  nn.Sequential(*list(network.disease_classifier.children()) , nn.Softmax(dim=1))
    def forward(self, input_data):
        celltype = self.celltype(self.feature(input_data))
        disease = self.disease(self.feature(input_data))
        return celltype, disease


def tune1(test_df, network, params):
    test_dat = Datasets(test_df)
    lr = params[0]
    n_epoch = params[1]
    batch_size = params[2]
    optimizer = torch.optim.Adam(network.parameters(), lr=lr)
    loss = nn.MSELoss()
    loss = loss.to(device)
    test_loader = DataLoader(dataset=test_dat,
                             batch_size=batch_size,
                             shuffle=True)
    for name, paras in network.encoder.named_parameters():
        paras.requires_grad = False
    network = network.to(device)
    for epoch in tqdm(range(n_epoch)):
        for batch in test_loader:
            expr = batch
            expr = expr.float()
            expr = expr.to(device)
            output = network(expr)
            err = loss(output, expr)
            optimizer.zero_grad()
            err.backward()
            optimizer.step()
    print('Finish Tuning1')
    return network


def tune2(test_df, network, params):
    test_dat = Datasets(test_df)
    lr = params[0]
    n_epoch = params[1]
    batch_size = params[2]
    optimizer = torch.optim.Adam(network.parameters(), lr=lr)
    loss = nn.MSELoss()
    loss = loss.to(device)
    test_loader = DataLoader(dataset=test_dat,
                             batch_size=batch_size,
                             shuffle=True)
    for name, paras in network.encoder.named_parameters():
        paras.requires_grad = True
    for name, paras in network.decoder.named_parameters():
        paras.requires_grad = False
    network = network.to(device)
    for epoch in tqdm(range(n_epoch)):
        for batch in test_loader:
            expr = batch
            expr = expr.float()
            expr = expr.to(device)
            output = network(expr)
            err = loss(output, expr)
            optimizer.zero_grad()
            err.backward()
            optimizer.step()
    print('Finish Tuning2')
    return network


def test(test_df, network, ct_dic, disease_dic ,disease):
    test_dat = Datasets(test_df)
    pred_prob = []
    ct_dic_rev = {v: k for k, v in ct_dic.items()}
    if disease:
        disease_dic_rev = {v: k for k, v in disease_dic.items()}
    test_loader = DataLoader(dataset=test_dat,
                             batch_size=test_df.shape[1],
                             shuffle=False)
    with torch.no_grad():
        pred_labels = []
        disease_labels = []
        for batch in test_loader:
            expr = batch
            expr = expr.float()
            expr = expr.to(device)
            if disease:
                class_output, disease_output = network(expr)
                disease_labels.append(
                disease_output.argmax(dim=1).cpu().numpy().tolist())
            else:
                class_output = network(expr)
            pred_labels.append(
                class_output.argmax(dim=1).cpu().numpy().tolist())
            pred_prob.append(class_output.cpu().numpy())
        pred_labels = [ct_dic_rev[i] for item in pred_labels for i in item]
        pred_prob = pd.DataFrame(reduce(pd.concat, pred_prob))
        pred_prob.index = test_df.columns
        pred_prob.columns = ct_dic.keys()
        if disease:
            disease_labels = [disease_dic_rev[i] for item in disease_labels for i in item]
            return pred_labels, pred_prob, disease_labels
        else:
            return pred_labels, pred_prob
    


def query_predict(query_expr, model, path_out, outprefix, disease):

    try:
        os.makedirs(path_out)
    except OSError:
        pass

    print('Loading data')
    query_expr = read_expr(query_expr)
    with open(model.replace('params', 'meta').replace('pt', 'pkl'), 'rb') as f:
        meta = pickle.load(f)
    if disease:
        genes, disease_dic, ct_dic = meta['genes'], meta['cellsources'], meta['celltypes']
    else:
        genes, ct_dic = meta['genes'], meta['celltypes']
    query_expr = merge(genes, query_expr)
    nfeatures, nct = len(genes), len(ct_dic)
    network = torch.load(model, map_location=device)
    if disease:
        network = Disease_Classifier(network).to(device)
        pred_labels, pred_prob, disease_labels = test(query_expr, network, ct_dic, disease_dic, disease)
        disease_labels = pd.DataFrame({'Cell':query_expr.columns,'Prediction':disease_labels})
        pd.DataFrame(disease_labels).to_csv(path_out + '/' + outprefix +
                                     '_cellsources.txt',
                                     index=False,
                                     header=True,
                                     sep='\t')
    else:
        network = Autoencoder(network, nfeatures, nct)
        print('Fine-tuning1')
        network = tune1(query_expr, network, params_tune1)
        print('Fine-tuning2')
        network = tune2(query_expr, network, params_tune2)
        network = Normal_Classifier(network).to(device) 
        pred_labels, pred_prob = test(query_expr, network, ct_dic, None, disease)

    pd.DataFrame(pred_labels).to_csv(path_out + '/' + outprefix +
                                     '_predictions.txt',
                                     index=False,
                                     header=False,
                                     sep='\t')
    pd.DataFrame(pred_prob).to_csv(path_out + '/' + outprefix +
                                   '_probability.txt',
                                   index=True,
                                   header=True,
                                   sep='\t')
    print('Finish Prediction')


def predict_downstream(seurat, prob_cutoff, path_out,
                       outprefix):
    print('Begin downstream analysis')
    cmd = 'Rscript ' + os.path.split(
        os.path.abspath(__file__)
    )[0] + '/downstream.R ' + ' --seurat ' + seurat + ' --pred_label ' + path_out + '/' + outprefix + '_predictions.txt' + ' --pred_prob ' + path_out + '/' + outprefix + '_probability.txt' + ' --prob_cutoff ' + str(
            prob_cutoff
        ) + ' --path_out ' + path_out + '  --outprefix ' + outprefix
    os.system(cmd)
    print('Finish downstream analysis')
