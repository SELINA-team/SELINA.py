# Run SELINA

## Pre-train

If you want to train a model with your own reference files, you can use the following command:

```
usage: selina train [-h] --path-in PATH_IN --path-out PATH_OUT
                    [--outprefix OUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --path-in PATH_IN     File path of training datasets.
  --path-out PATH_OUT   File path of the output model.
  --outprefix OUTPREFIX Prefix of the output files. DEFAULT: pre-trained
```

In this step, two output files used in the next step will be generated.

- `pre-trained_params.pt` : a file containing all parameters of the pre-trained model
- `pre-trained_meta.pkl` : a file containing the cell types and genes of the reference data

## Predict

Here you can choose one model from our pre-trained models (available on [SELINA models](https://github.com/wanglabtongji/SELINA_reference)) or the model trained by yourself to annotate the query data.

```
usage: selina predict [-h] --mode {single,cluster} --query-expr QUERY_EXPR
                      --model MODEL --seurat SEURAT [--disease]
                      [--cell-cutoff CELL_CUTOFF] [--prob-cutoff PROB_CUTOFF]
                      [--path-out PATH_OUT] [--outprefix OUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit

Arguments for input.:
  --mode {single,cluster}
                        Single-cell level input or cluster level input.
  --query-expr QUERY_EXPR
                        File path of query data matrix.
  --model MODEL         File path of the pre-trained model.
  --seurat SEURAT       Path of seurat object.
  --disease             The input data is from disease condition.

Cutoff for downstream analysis:
  --cell-cutoff CELL_CUTOFF
                        Cutoff for cells with the same cell type in 10 nearest
                        neighbor cells(only used when the input is sinle-cell
                        level expression matrix). DEFAULT: 5.
  --prob-cutoff PROB_CUTOFF
                        Cutoff for prediction probability. DEFAULT: 0.9.

Output Arguments:
  --path-out PATH_OUT   File path of the output files.
  --outprefix OUTPREFIX
                        Prefix of the output files. DEFAULT: query
```

This step will output four files:

- `query_predictions.txt`: predicted cell type for each cell in the query data

- `query_probability.txt`: probability of cells predicted as each of the reference cell types

- `query_pred.png`: UMAP plot with cell type labels

- `query_DiffGenes.tsv`: differentially expressed genes for each cell type

if the input data is one single-cell level expression matrix, four additional files will be generated, which are:

- `query_cluster_prob.png`: box plot indicating the prediction probability distribution of cells in each cluster

- `query_prob.txt`: text file containing the prediction probability of each cell

- `query_unknown_percent.png`: bar plot indicating the percentage of cells that are assigned unknown in each cluster

- `query_unknown_percent.txt`: text file representing the percentage of cells that are assigned unknown in each cluster
