# Run SELINA

## Pre-train

You can follow the usage instruction listed below if you want to train a model with your own reference files. Note that you should include `--disease` in your command if you want to train models for disease data.

```
usage: selina train [-h] --path-in PATH_IN --path-out PATH_OUT
                    [--outprefix OUTPREFIX] [--disease]

optional arguments:
  -h, --help            show this help message and exit
  --path-in PATH_IN     File path of training datasets.
  --path-out PATH_OUT   File path of the output model.
  --outprefix OUTPREFIX
                        Prefix of the output files. DEFAULT: pre-trained
  --disease             This flag should be used when the data is in
                        disease condition
```

In this step, two output files used in the next step will be generated.

- `pre-trained_params.pt` : a file containing all parameters of the pre-trained model
- `pre-trained_meta.pkl` : a file containing the cell types and genes of the reference data

## Predict

Here you can choose one model from our pre-trained models (available on [SELINA models](https://github.com/SELINA-team/SELINA-reference)) or the model trained by yourself to annotate the query data. These pre-trained models are divided into two categories, of which one is for normal dataset prediction, and another one is for disease datasets prediction. Currently the disease models only cover the non-small-cell lung carcinoma, type 2 diabetes and Alzheimer's disease, which were used to evaluate the performance of SELINA in our paper.

Since the expression profiles of disease data may be more complicated than normal data, we removed the fine-tuning step when predicting for the disease data, which can be achieved by adding `--disease` to the command.

```
usage: selina predict [-h] --query-expr QUERY_EXPR --model MODEL --seurat
                      SEURAT [--disease] [--prob-cutoff PROB_CUTOFF]
                      --path-out PATH_OUT [--outprefix OUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit

Arguments for input:
  --query-expr QUERY_EXPR
                        File path of the query data matrix.
  --model MODEL         File path of the pre-trained model.
  --seurat SEURAT       File path of the seurat object.
  --disease             This flag should be used when the data is in some
                        disease condition

Cutoff for downstream analysis:
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
