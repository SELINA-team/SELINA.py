# Run SELINA

## Train

If you want to train a model with your own reference files, you can use the following command:

```
usage: selina train [-h] --path_in PATH_IN --path_out PATH_OUT
                     [--outprefix OUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --path_in PATH_IN     File path of training datasets.
  --path_out PATH_OUT   File path of the output model.
  --outprefix OUTPREFIX Prefix of the output files. DEFAULT: pre-trained
```

In this step, two output files used in the next step will be generated.

- `pre-trained_params.pt` : a file containing all parameters of the trained model
- `pre-trained_meta.pkl` : a file containing the cell types and genes of the reference data

## Predict

Here you can choose to use our pre-trained models (available on [SELINA models](https://github.com/wanglabtongji/SELINA_reference)) or the model trained by yourself to annotate the query data.

```
usage: selina predict [-h] --mode {single,cluster} --input INPUT --model MODEL
                       --path_out PATH_OUT [--outprefix OUTPREFIX] --plot
                       {True,False} [--rds RDS]

optional arguments:
  -h, --help            show this help message and exit

Arguments for prediction.:
  --mode {single,cluster}
                        Single-cell level input or cluster level input.
  --input INPUT         File path of query data.
  --model MODEL         File path of the pre-trained model.
  --path_out PATH_OUT   File path of the output files.
  --outprefix OUTPREFIX
                        Prefix of the output files. DEFAULT: query

Arguments for plot:
  --plot {True,False}   Whether to generate umap plot. DEFAULT: False
  --rds RDS             File path of rds file.
```

This step will output three files:

- `query_predictions.txt` : predicted cell type for each cell in the query data(choose the cell type corresponding to the max probablity as the default prediction results)
- `query_probability.txt` : probablity of cells predicted as each of the reference cell types
- `query_pred.png` : a umap png file with cell type annotation on it
