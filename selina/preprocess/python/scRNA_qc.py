import os
import numpy as np
import scipy.sparse as sp_sparse
from selina.preprocess.python.scRNA_data import *

def Filter(rawmatrix, feature, barcode, count_cutoff, gene_cutoff, cell_cutoff,
           outprefix, assembly):
    count_per_cell = np.asarray(rawmatrix.sum(axis=0))
    genes_per_cell = np.asarray((rawmatrix > 0).sum(axis=0))
    count_gene = np.concatenate((count_per_cell, genes_per_cell), axis=0)
    statfile = outprefix + "_count_gene_stat.txt"
    with open(statfile, "w") as stat_out:
        header = "Cell\tCount\tGene\n"
        stat_out.write(header)
        for i in range(count_gene.shape[1]):
            stat_list = count_gene[0:2, i].tolist()
            stat_list = [str(int(j)) for j in stat_list]
            stat_out.write(barcode[i] + "\t" + "\t".join(stat_list) + "\n")

    passed_cell = np.logical_and(count_per_cell > count_cutoff,
                                 genes_per_cell > gene_cutoff)

    cells_per_gene = np.asarray((rawmatrix > 0).sum(axis=1))
    passed_gene = cells_per_gene > cell_cutoff
    passed_gene = np.transpose(passed_gene)

    passed_cell_matrix = rawmatrix[np.where(passed_gene.flatten())[0], :]
    passed_cell_matrix = passed_cell_matrix[:,
                                            np.where(passed_cell.flatten())[0]]

    passed_barcodes = np.array(barcode)[passed_cell.tolist()[0]].tolist()
    passed_genes = np.array(feature)[passed_gene.tolist()[0]].tolist()

    write_10X_h5(outprefix + "_filtered_gene_count.h5",
                 matrix=passed_cell_matrix,
                 features=passed_genes,
                 barcodes=passed_barcodes,
                 genome=assembly,
                 datatype='Gene')


def scrna_qc(directory, outprefix, fileformat, matrix, separator, feature,
             gene_column, barcode, count_cutoff, gene_cutoff, cell_cutoff,
             assembly):

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    if fileformat == "plain":
        matrix_dict = read_count(matrix, separator)
        rawmatrix = matrix_dict["matrix"]
        rawmatrix = sp_sparse.csc_matrix(rawmatrix, dtype=numpy.float32)
        features = matrix_dict["features"]
        barcodes = matrix_dict["barcodes"]

        if features[0][0] == '"' or features[0][0] == "'":
            features = [i[1:(len(i) - 1)] for i in features]
        if barcodes[0][0] == '"' or barcodes[0][0] == "'":
            barcodes = [i[1:(len(i) - 1)] for i in barcodes]

        features = [i.split('_')[-1] for i in features]
        barcodes = [
            i.split('/')[-1].split('.genes.results')[0] for i in barcodes
        ]

    elif fileformat == "h5":
        scrna_count = read_10X_h5(matrix)
        rawmatrix = scrna_count.matrix
        features = scrna_count.names.tolist()
        barcodes = scrna_count.barcodes.tolist()

        if type(features[0]) == bytes:
            features = [i.decode() for i in features]
        if type(barcodes[0]) == bytes:
            barcodes = [i.decode() for i in barcodes]

    elif fileformat == "mtx":
        matrix_dict = read_10X_mtx(matrix_file=matrix,
                                   feature_file=feature,
                                   barcode_file=barcode,
                                   datatype="Gene",
                                   gene_column=gene_column)
        rawmatrix = matrix_dict["matrix"]
        features = matrix_dict["features"]
        barcodes = matrix_dict["barcodes"]

    filename = os.path.join(directory, outprefix)

    Filter(rawmatrix=rawmatrix,
            feature=features,
            barcode=barcodes,
            count_cutoff=count_cutoff,
            gene_cutoff=gene_cutoff,
            cell_cutoff=cell_cutoff,
            outprefix=filename,
            assembly=assembly)
