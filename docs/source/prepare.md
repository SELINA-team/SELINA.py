# Prepare

## Installation

SELINA is available for macOS, Linux and Windows.
To use SELINA, you should first build a conda environment.

```
conda create -n Selina
conda activate Selina
```

All the dependency packages will be installed simultaneously with the following commands except for the presto which can be found at [presto](https://github.com/immunogenomics/presto). The R package devtools used in the installation of presto has been included in SELINA, so you do not need to install devtools once again.

```
conda install -c conda-forge -c r -c bioconda -c pfren selina
```

Note that if you have gpu on your device and want to use it, you should additionally run the following command to install cudatoolkit and the paired version of pytorch and cudnn.

```
conda install pytorch cudatoolkit cudnn
```

## Prepare data

### Preprocess of training data

Before you start to run SELINA, make sure you have prepared the reference datasets with the format shown as below. For each reference dataset, you should have 2 paired files, one is named as xx_expr.txt which contains the gene expression profile, and the other is named as xx_meta.txt containing the meta information for each cell. For the expression profile, the first column is gene name which is followed by expression of each cell, and the gene names should be symbol names and matched with hg38. For the meta file, the first column is cell type of each cell, and the second column is the sequencing platform of each cell. If you want to train data with normal cells and disease cells mixed, you should additionally add one column named Disease which presents the cell source(Normal/Name of any disease). Note that if you choose to use our pre-trained models, this step can be skipped. We recommend that you use the original count data as the reference data.

#### File tree

```
reference/
├── train1_expr.txt
├── train1_meta.txt
├── train2_expr.txt
└── train2_meta.txt
```

#### Expression matrix

| Gene  | Cell1 | Cell2 |
| :---: | :---: | :---: |
| NOC2L |   7   |   3   |
| ISG15 |  10   |   2   |

#### Meta files

##### Meta file for normal data

| Celltype | Platform |
| :------: | :------: |
|   CD8T   |   10x    |
|   CD4T   |   10x    |

##### Meta file for disease data

| Celltype | Platform | Disease |
| :------: | :------: | :-----: |
|   CD8T   |   10x    |   T2D   |
|   CD4T   |   10x    | Normal  |

### Preprocess of query data

This step is to normalize, convert the genes to version hg38 and symbol names, perform dimension reduction and clustering for your data. SELINA supports 3 formats of input: `plain`,`h5` and `mtx`. The gene by cell matrix is in plain format. The full list of preprocessing commands is shown as below:

```
usage: selina preprocess [-h] --format {h5,mtx,plain} [--matrix MATRIX]
                         [--separator {tab,space,comma}] [--feature FEATURE]
                         [--gene-column GENE_COLUMN] [--barcode BARCODE]
                         [--gene-idtype {symbol,ensembl}]
                         [--assembly {GRCh38,GRCh37}]
                         [--count-cutoff COUNT_CUTOFF]
                         [--gene-cutoff GENE_CUTOFF]
                         [--cell-cutoff CELL_CUTOFF] [--mito]
                         [--mito-cutoff MITO_CUTOFF]
                         [--variable-genes VARIABLE_GENES] [--npcs NPCS]
                         [--cluster-res CLUSTER_RES] --directory DIRECTORY
                         [--outprefix OUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit

Input files arguments:
  --format {h5,mtx,plain}
                        Format of the count matrix file.
  --matrix MATRIX       Location of count matrix file. If the format is 'h5'
                        or 'plain', users need to specify the name of the
                        count matrix file.If the format is 'mtx', the
                        'matrix' should be the name of .mtx formatted matrix
                        file, such as 'matrix.mtx'.
  --separator {tab,space,comma}
                        The separating character (only for the format of
                        'plain').Values on each line of the plain matrix file
                        will be separated by the character. DEFAULT: tab.
  --feature FEATURE     Location of feature file (required for the format of
                        'mtx'). Features correspond to row indices of count
                        matrix. DEFAULT: features.tsv.
  --gene-column GENE_COLUMN
                        If the format is 'mtx', please specify which column
                        of the feature file to use for gene names. DEFAULT:
                        2.
  --barcode BARCODE     Location of barcode file (required for the format of
                        'mtx'). Cell barcodes correspond to column indices of
                        count matrix. DEFAULT: barcodes.tsv.
  --gene-idtype {symbol,ensembl}
                        Type of gene name, 'symbol' for gene symbol and
                        'ensembl' for ensembl id. DEFAULT: symbol.
  --assembly {GRCh38,GRCh37}
                        Assembly (GRCh38/hg38 and GRCh37/hg19). DEFAULT:
                        GRCh38.

Quality control arguments:
  --count-cutoff COUNT_CUTOFF
                        Cutoff for the number of count in each cell. DEFAULT:
                        1000.
  --gene-cutoff GENE_CUTOFF
                        Cutoff for the number of genes included in each cell.
                        DEFAULT: 500.
  --cell-cutoff CELL_CUTOFF
                        Cutoff for the number of cells covered by each gene.
                        DEFAULT: 10.
  --mito                This flag should be used when you want to filter out
                        cells with a high percentage of mitochondria genes.
  --mito-cutoff MITO_CUTOFF
                        Cutoff for the percentage of mitochondria genes in
                        each cell. DEFAULT: 0.2.

Process arguments:
  --variable-genes VARIABLE_GENES
                        Number of variable genes used in PCA. DEFAULT: 2000.
  --npcs NPCS           Number of dimensions after PCA. DEFAULT: 30.
  --cluster-res CLUSTER_RES
                        Clustering resolution. DEFAULT: 0.6.

Output arguments:
  --directory DIRECTORY
                        Path to the directory where the result file shall be
                        stored.
  --outprefix OUTPREFIX
                        Prefix of output files. DEFAULT: query.
```

In this step four output files will be generated:

- `query_res.rds`: a seurat object storing the normalized data, dimension reduction and clustering results

- `query_expr.txt`: expression matrix of query data with the first column as genes

- `query_cluster.png`: UMAP plot with cluster labels

- `query_cluster_DiffGenes.tsv`: differentially expressed genes for each cluster
