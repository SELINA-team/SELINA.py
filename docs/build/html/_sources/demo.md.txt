# Demo

This is an example to show you how to run selina step by step.

## 1. Data

The data used in this vignette are orgnized as the following directory tree shows.

```
├── query_data
│   └── GSE139534_Lung-Adult_10462_gene_count.h5
├── reference_data
│   ├── GSE123405_Lung-Adult_7786_expr.txt
│   ├── GSE123405_Lung-Adult_7786_meta.txt
│   ├── GSE130148_Lung-Adult_2306_expr.txt
│   ├── GSE130148_Lung-Adult_2306_meta.txt
│   ├── GSE134355_Lung-Adult_5672_expr.txt
│   ├── GSE134355_Lung-Adult_5672_meta.txt
│   ├── GSE134355_Lung-Adult_8426_expr.txt
│   ├── GSE134355_Lung-Adult_8426_meta.txt
│   ├── GSE134355_Lung-Adult_9113_expr.txt
│   ├── GSE134355_Lung-Adult_9113_meta.txt
│   ├── GSE146981_Lung-Adult_27015_expr.txt
│   └── GSE146981_Lung-Adult_27015_meta.txt
└── res
```

The query data used here is from [mTORC1 activation in lung mesenchyme drives sex- and age-dependent pulmonary structure and function decline](https://www.nature.com/articles/s41467-020-18979-4). The citation papers of the reference data are listed in the following table.

|          Dataset           |                                                                                                        Study                                                                                                        |
| :------------------------: | :-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| GSE123405_Lung-Adult_7786  |               [Cell-specific expression of lung disease risk-related genes in the human small airway epithelium](https://respiratory-research.biomedcentral.com/articles/10.1186/s12931-020-01442-9)                |
| GSE130148_Lung-Adult_2306  |                                     [A cellular census of human lungs identifies novel cell states in health and in asthma](https://www.nature.com/articles/s41591-019-0468-5)                                      |
| GSE134355_Lung-Adult_5672  |                                                  [Construction of a human cell landscape at single-cell level](https://www.nature.com/articles/s41586-020-2157-4)                                                   |
| GSE134355_Lung-Adult_8426  |                                                  [Construction of a human cell landscape at single-cell level](https://www.nature.com/articles/s41586-020-2157-4)                                                   |
| GSE134355_Lung-Adult_9113  |                                                  [Construction of a human cell landscape at single-cell level](https://www.nature.com/articles/s41586-020-2157-4)                                                   |
| GSE146981_Lung-Adult_27015 | [Senescence of Alveolar Type 2 Cells Drives Progressive Pulmonary Fibrosis](https://www.atsjournals.org/doi/10.1164/rccm.202004-1274OC?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) |

Specifically, the reference datasets are stored in the reference_data folder, and for each dataset there are two corresponding files used in the pretraining stage; one with `expr` provides the expression profile of each dataset, and the file with `meta` presents the cell type and sequencing platform of each cell. Here we take GSE123405_Lung-Adult_7786 as an example to show you the detail information.

```
GSE123405_Lung-Adult_7786_expr.txt
```

```
Gene GSM3502715@GTTAGCTTAATT GSM3502715@TCCCCTCTTCGC GSM3502715@TATTGTTTTACN GSM3502715@GACTTATTTATA
A2M 0 0 0 0
A4GALT 1 2 1 0
AAAS 0 0 3 0
AACS 0 2 0 0
AADAC 0 0 0 0
```

```
GSE123405_Lung-Adult_7786_meta.txt
```

```
Celltype Platform
Club Drop-seq
Club Drop-seq
Club Drop-seq
Ciliated Drop-seq
```

## 2. Preprocess of query data

The following command is to normalize, match the assembly version with the reference data, perform dimension reduction for your data. We support 3 formats of input: `plain`,`h5` and `mtx`. The plain format is a gene by cell matrix.

```
selina preprocess --format h5  --matrix query_data/GSE139534_Lung-Adult_10462_gene_count.h5 --gene-idtype symbol --assembly GRCh38 --count-cutoff 1000 --gene-cutoff 500 --cell-cutoff 10 --directory res --outprefix query --mode single
```

```
Normalization and identify variable genes ...
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Regressing out nCount_RNA
  |======================================================================| 100%
Centering and scaling data matrix
  |======================================================================| 100%
PCA analysis ...
PC_ 1
Positive:  TMSB4X, SRGN, VIM, PECAM1, HLA-E, EMP3, BST2, GMFG, HLA-C, GRN
	   TXNIP, TMSB10, GSTO1, SH3BGRL3, ENG, HLA-A, HLA-B, IFI27, GPX3, KLF4
	   FKBP1A, AKR1B1, GYPC, TIMP1, PLIN2, GIMAP4, IFNGR1, CD68, THBD, ICAM2
Negative:  C20orf85, FAM183A, C1orf194, RSPH1, TMEM190, C9orf24, LRRIQ1, C11orf88, SNTN, C9orf116
	   DNAAF1, FAM92B, CCDC78, DYNLRB2, PIFO, CAPS, MS4A8, C5orf49, CAPSL, CFAP126
	   ZMYND10, TPPP3, AC013264.1, MORN5, TSPAN1, RSPH9, DYDC2, CFAP53, CCDC146, DRC3
PC_ 2
Positive:  SERPINA1, MS4A7, LGALS3BP, VSIG4, FTL, MARCO, SPI1, FBP1, ACP5, C1orf162
	   CD163, CAPG, ALOX5AP, MCEMP1, OLR1, AIF1, EVI2B, CD37, GPNMB, CD53
	   TREM1, CSTA, FCGR3A, SNX10, MS4A4A, APOC1, C1QC, CD68, CFD, MSR1
Negative:  TIMP3, CLDN5, SPARCL1, RAMP2, PCAT19, SRPX, VWF, PRSS23, IGFBP7, HYAL2
	   MT1M, NPDC1, CLEC14A, PTPRB, CAV1, GIMAP7, SLCO2A1, ID1, CRIP2, FAM107A
	   GNG11, SLC9A3R2, TM4SF1, CAVIN1, CALCRL, TSPAN7, EPAS1, TINAGL1, MT1X, A2M
PC_ 3
Positive:  SRGN, B2M, VIM, CYBA, GMFG, MS4A7, CD163, LGALS1, LAPTM5, VSIG4
	   C1QC, MARCO, MRC1, SPI1, ALOX5AP, C1orf162, SLC11A1, TYROBP, INHBA, AIF1
	   MSR1, OLR1, PDK4, GPNMB, ALOX5, CFD, SNX10, IFI27, CD68, PECAM1
Negative:  SFTA2, HOPX, PEBP4, SFTPD, NAPSA, S100A14, SFTA3, PGC, MUC1, SFTPB
	   CLDN18, SDR16C5, CXCL17, SOD3, SLC34A2, GKN2, SLC39A8, LAMP3, SFTPA2, C16orf89
	   ABCA3, SFTPA1, SELENOP, TMEM97, SFN, PLA2G1B, WIF1, LPCAT1, MFSD2A, FASN
PC_ 4
Positive:  CRYAB, ACTA2, WISP2, PLN, PPP1R14A, ACTG2, MYL9, NEXN, TPM2, MYH11
	   CNN1, LMOD1, ADH1B, DES, COL6A2, BGN, MFAP4, PDGFRB, PRELP, NTM
	   TAGLN, NDUFA4L2, KCNMB1, FXYD1, NOTCH3, MAP1B, TNNT2, FMO2, SCEL, SPOCK2
Negative:  CD74, NPC2, HLA-DRB1, WIF1, LAMP3, SFTPA1, SFTPA2, MFSD2A, CTSH, LRRK2
	   SERPINA1, TCIM, PLA2G1B, NEAT1, SOCS3, SFTPC, TMEM97, FABP5, HLA-DRA, HLA-DPA1
	   CXCL17, FOS, LIFR, AQP1, FCN3, JUNB, CYB5A, HMGCS1, FASN, ABCA3
PC_ 5
Positive:  NTM, AGER, SCEL, RTKN2, CEACAM6, SPOCK2, MS4A15, ITLN2, CLIC5, GGTLC1
	   GGT1, ABCA7, MYRF, SEMA3B, TNNC1, ANKRD29, KLK11, FMO2, LAMA3, COL12A1
	   C19orf33, KRT7, PLCXD2, SLC1A1, GPRC5A, TACSTD2, CST6, ARHGEF26, SCNN1B, IL32
Negative:  WISP2, ACTA2, PLN, ACTG2, MYH11, LMOD1, CNN1, ADH1B, PPP1R14A, C11orf96
	   DES, COL6A2, BGN, TPM2, PDGFRB, IGFBP6, MFAP4, WIF1, PRELP, FXYD1
	   LAMP3, SOD3, NDUFA4L2, MFSD2A, LRRK2, SFTPA1, MAP1B, TNNT2, PLA2G1B, SFTPA2
UMAP analysis ...
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
10:59:39 UMAP embedding parameters a = 0.9922 b = 1.112
10:59:39 Read 10453 rows and found 28 numeric columns
10:59:39 Using Annoy for neighbor search, n_neighbors = 30
10:59:39 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:59:40 Writing NN index file to temp file /tmp/Rtmp0WAzsB/file2d504385ba56
10:59:40 Searching Annoy index using 1 thread, search_k = 3000
10:59:43 Annoy recall = 100%
10:59:43 Commencing smooth kNN distance calibration using 1 thread
10:59:43 Initializing from normalized Laplacian + noise
10:59:44 Commencing optimization for 200 epochs, with 457718 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:59:48 Optimization finished
Computing nearest neighbor graph
Computing SNN
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 10453
Number of edges: 419653

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8941
Number of communities: 24
Elapsed time: 1 seconds
```

This step will output two files in the res folder.

- `1. query_res.rds` : a seurat object with gene expression profile and dimension reduction result

```
An object of class Seurat
17547 features across 10453 samples within 1 assay
Active assay: RNA (17547 features, 2000 variable features)
2 dimensional reductions calculated: pca, umap
```

- `2. query_single_expr.txt` : expression matrix of query data for the prediction step

```
Gene GSM4143262@AAACCTGCATCAGTCA-1 GSM4143262@AAACCTGTCGGGAGTA-1 GSM4143262@AAAGTAGTCGGAATCT-1 GSM4143262@AACCATGCACGACGAA-1
<chr> <dbl> <dbl> <dbl> <dbl>
AL669831.5 0 0 0 0
FAM87B 0 0 0 0
LINC00115 0 0 0 0
FAM41C 0 0 0 0
AL645608.1 0 0 0 0
```

## 3. Train

This step is to train a model using the reference data listed above.

```
selina train --path_in reference_data --path_out res --outprefix pre-trained
```

```
Loading data
100% |██████████████████████████████████████████████████| Reading data [done]
100% |██████████████████████████████████████████████████| Reading data [done]
100% |██████████████████████████████████████████████████| Reading data [done]
100% |██████████████████████████████████████████████████| Reading data [done]
100% |██████████████████████████████████████████████████| Reading data [done]
100% |██████████████████████████████████████████████████| Reading data [done]
Begin training
100%|████████████████████████████████████████████████████████████████████████| 50/50 [11:59<00:00, 14.39s/it]
Finish Training
All done
```

In this step, two output files used in the next step will be generated.

- `1. pre-trained_params.pt` : a file containing all parameters of the trained model

- `2. pre-trained_meta.pkl` : a file containing the cell types and genes of the reference data

```
with open('res/pre-trained_meta.pkl','rb') as f:
meta = pickle.load(f)
meta['genes'][1:5]
array(['A4GALT', 'AAAS', 'AACS', 'AADAC'], dtype='<U12')

meta['celltypes'].keys()
dict_keys(['Mucous', 'CD8T', 'Fibroblast', 'Macrophage', 'Chondrocyte', 'DC', 'AT2', 'Endothelial_lv3', 'NK', 'Neutrophil', 'Alveolar Bipotent', 'AT1', 'Lymphatic Endothelial', 'Muscle', 'Mast', 'Secretory Epithelial', 'Club', 'Basal', 'Alveolar Bipotent Progenitor', 'Goblet', 'CD4T', 'Megakaryocyte', 'Intermediate(Club-Basal)', 'Ionocyte', 'B', 'Neuroendocrine', 'Ciliated', 'Monocyte'])
```

## 4. Predict

Here you can choose to use our pretrained models (available on [SELINA models](https://github.com/wanglabtongji/SELINA_reference)) or the model trained by yourself to annotate the query data.

```
selina predict --mode single --input res/query_single_expr.txt --model res/pre-trained_params.pt --path_out res --outprefix query --plot True --rds res/query_res.rds
```

```
Loading data
100% |██████████████████████████████████████████████████| Reading data [done]
Fine-tuning1
100%|████████████████████████████████████████████████████████████████████████| 50/50 [01:26<00:00, 1.73s/it]
Finish Tuning1
Fine-tuning2
100%|████████████████████████████████████████████████████████████████████████| 10/10 [00:24<00:00, 2.49s/it]
Finish Tuning2
Finish Prediction
Plotting
Finish plotting
```

This step will output three files:

- `1. query_predictions.txt` : predicted cell type for each cell in the query data(choose the cell type corresponding to the max probablity as the default prediction results)

```
Cell Prediction
GSM4143262@AAACCTGCATCAGTCA-1 AT1
GSM4143262@AAACCTGTCGGGAGTA-1 AT1
GSM4143262@AAAGTAGTCGGAATCT-1 AT1
GSM4143262@AACCATGCACGACGAA-1 AT1
GSM4143262@AACGTTGAGGCGCTCT-1 AT1
```

- `2. query_probability.txt` : probablity of cells predicted as each of the reference cell types

```
Mucous CD8T Fibroblast Macrophage
GSM4143262@AAACCTGCATCAGTCA-1 9.528222e-13 3.493414e-26 8.598993e-16 2.151410e-10
GSM4143262@AAACCTGTCGGGAGTA-1 6.421099e-15 1.029570e-29 4.329293e-17 4.808651e-14
GSM4143262@AAAGTAGTCGGAATCT-1 2.036562e-18 2.974594e-34 3.507774e-21 3.168965e-19
GSM4143262@AACCATGCACGACGAA-1 6.580991e-16 4.516618e-29 9.358461e-18 5.502941e-14
GSM4143262@AACGTTGAGGCGCTCT-1 2.829751e-13 5.663002e-25 3.487227e-16 1.142975e-12
```

- `3. query_pred.png` : a umap png file with cell type annotation on it

![image](./_images/query_pred.png)
