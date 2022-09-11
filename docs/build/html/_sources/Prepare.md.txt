# Prepare

## Installation

SELINA is supported for macOS, Linux and Windows.
To use SELINA, you should first build a conda environment.

```
conda create -n Selina
conda activate Selina
```

Then you can install SELINA using the following command(All the dependency packages will be installed simultaneously).

```
conda install -c pfren selina -c conda-forge -c r
```

Note that if you have gpu on your device, you should additionally run the following command after the above commands are executed.

```
conda install pytorch cudatoolkit
```

## Prepare data

### Preprocess of training data

Before you start to run SELINA, make sure you have prepared the proper reference datasets, the format is shown as below. For each reference datasets, you should have 2 paired files, one is named as xx_expr.txt, this file contains the gene expression, and the other is named as xx_meta.txt. For the expression profile, the first column is gene which is followed by expression of each cell.
For the meta file, the first column is celltype of each cell, and the second column is platform of each cell. Note that if you choose to use our pretarined models, this step can be skipped.

```
reference/
├── train1_expr.txt
├── train1_meta.txt
├── train2_expr.txt
└── train2_meta.txt
```

| Gene  | Cell1 | Cell2 |
| :---: | :---: | :---: |
| NOC2L |   7   |   3   |
| ISG15 |  10   |   2   |

| Celltype | Platform |
| :------: | :------: |
|   CD8T   |   10x    |
|   CD4T   |   10x    |

### Preprocess of query data

In addition to the training data, you also need to preprocess the query data. This step is to normalize, match the assembly version with the reference data, perform dimension reduction for your data. We support 3 formats of input: `plain`,`h5` and `mtx`. The plain format is a gene by cell matrix. The full list of preprocessing commands is shown as below:

```
usage: selina preprocess [-h] [--format {h5,mtx,plain}] [--matrix MATRIX]
                          [--separator {tab,space,comma}] [--feature FEATURE]
                          [--gene-column GENE_COLUMN]
                          [--gene-idtype {symbol,ensembl}] [--barcode BARCODE]
                          [--assembly {GRCh38,GRCh37}]
                          [--count-cutoff COUNT_CUTOFF]
                          [--gene-cutoff GENE_CUTOFF]
                          [--cell-cutoff CELL_CUTOFF] [--directory DIRECTORY]
                          [--outprefix OUTPREFIX] --mode {single,cluster,both}

optional arguments:
  -h, --help            show this help message and exit

Input files arguments:
  --format {h5,mtx,plain}
                        Format of the count matrix file.
  --matrix MATRIX       Location of count matrix file. If the format is 'h5'
                        or 'plain', users need to specify the name of the
                        count matrix file.If the format is 'mtx', the 'matrix'
                        should be the name of .mtx formatted matrix file, such
                        as 'matrix.mtx'.
  --separator {tab,space,comma}
                        The separating character (only for the format of
                        'plain').Values on each line of the plain matrix file
                        will be separated by the character. DEFAULT: tab.
  --feature FEATURE     Location of feature file (required for the format of
                        'mtx'). Features correspond to row indices of count
                        matrix. DEFAULT: features.tsv.
  --gene-column GENE_COLUMN
                        If the format is 'mtx', please specify which column of
                        the feature file to use for gene names. DEFAULT: 2.
  --gene-idtype {symbol,ensembl}
                        Type of gene name, 'symbol' for gene symbol and
                        'ensembl' for ensembl id. DEFAULT: symbol.
  --barcode BARCODE     Location of barcode file (required for the format of
                        'mtx'). Cell barcodes correspond to column indices of
                        count matrix. DEFAULT: barcodes.tsv.
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

Output arguments:
  --directory DIRECTORY
                        Path to the directory where the result file shall be
                        stored. DEFAULT: preprocess.
  --outprefix OUTPREFIX
                        Prefix of output files. DEFAULT: query.
  --mode {single,cluster,both}
                        Output expression file for prediction. single: single-
                        cell level. cluster: cluster level. both: output both
                        the single-cell level and cluster level expression
                        profiles
```

Note that you must choose the mode for the returned expression profiles. In this step two output files will be generated:

- `query_res.rds` : a seurat object with gene expression profile and dimension reduction result
- `query_{single/cluster}_expr.txt` : expression matrix of query data for the prediction step
