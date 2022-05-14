import numpy
import gzip
import h5py
import tables
import scipy.io
import collections
import scipy.sparse as sp_sparse


def read_count(count_file, separator="tab"):
    """Read count table as matrix."""

    if separator == "tab":
        sep = "\t"
    elif separator == "space":
        sep = " "
    elif separator == "comma":
        sep = ","
    else:
        raise Exception("Invalid separator!")

    infile = open(count_file, 'r').readlines()
    barcodes = infile[0].strip().split(sep)
    features = []
    matrix = []
    for line in infile[1:]:
        line = line.strip().split(sep)
        features.append(line[0])
        matrix.append([float(t) for t in line[1:]])
    if len(barcodes) == len(matrix[0]) + 1:
        barcodes = barcodes[1:]

    return {"matrix": matrix, "features": features, "barcodes": barcodes}

FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['ids','names','barcodes','matrix'])

def read_10X_h5(filename):
    """Read 10X HDF5 files, support both gene expression and peaks."""
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, 'matrix')
        except tables.NoSuchNodeError:
            print("Matrix group does not exist in this file.")
            return None
        feature_group = getattr(group, 'features')
        ids = getattr(feature_group, 'id').read()
        names = getattr(feature_group, 'name').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return FeatureBCMatrix(ids, names, barcodes, matrix)


def read_10X_mtx(matrix_file,
                 feature_file,
                 barcode_file,
                 datatype,
                 gene_column=2):
    """Convert 10x mtx as matrix."""

    matrix = scipy.io.mmread(matrix_file)
    matrix = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)

    if feature_file.split('.')[-1] == 'gz' or feature_file.split(
            '.')[-1] == 'gzip':
        feature_in = gzip.open(feature_file, "r")
    else:
        feature_in = open(feature_file, "r")
    features = feature_in.readlines()
    if datatype == "Peak":
        features = [
            "_".join(feature.strip().split("\t")[0:3]) for feature in features
        ]
    else:
        if type(features[0]) == str:
            features = [
                feature.strip().split("\t")[gene_column - 1]
                for feature in features
            ]
        if type(features[0]) == bytes:
            features = [
                feature.decode().strip().split("\t")[gene_column - 1]
                for feature in features
            ]

    if barcode_file.split('.')[-1] == 'gz' or barcode_file.split(
            '.')[-1] == 'gzip':
        barcode_in = gzip.open(barcode_file, "r")
    else:
        barcode_in = open(barcode_file, "r")
    barcodes = barcode_in.readlines()
    if type(barcodes[0]) == str:
        barcodes = [barcode.strip().split("\t")[0] for barcode in barcodes]
    if type(barcodes[0]) == bytes:
        barcodes = [
            barcode.decode().strip().split("\t")[0] for barcode in barcodes
        ]

    return {"matrix": matrix, "features": features, "barcodes": barcodes}


def write_10X_h5(filename,
                matrix,
                features,
                barcodes,
                genome='GRCh38',
                datatype='Peak'):
    """Write 10X HDF5 files, support both gene expression and peaks."""
    f = h5py.File(filename, 'w')
    if datatype == 'Peak':
        M = sp_sparse.csc_matrix(matrix, dtype=numpy.int8)
    else:
        M = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)
    B = numpy.array(barcodes, dtype='|S200')
    P = numpy.array(features, dtype='|S100')
    GM = numpy.array([genome] * len(features), dtype='|S10')
    FT = numpy.array([datatype] * len(features), dtype='|S100')
    AT = numpy.array(['genome'], dtype='|S10')
    mat = f.create_group('matrix')
    mat.create_dataset('barcodes', data=B)
    mat.create_dataset('data', data=M.data)
    mat.create_dataset('indices', data=M.indices)
    mat.create_dataset('indptr', data=M.indptr)
    mat.create_dataset('shape', data=M.shape)
    fet = mat.create_group('features')
    fet.create_dataset('_all_tag_keys', data=AT)
    fet.create_dataset('feature_type', data=FT)
    fet.create_dataset('genome', data=GM)
    fet.create_dataset('id', data=P)
    fet.create_dataset('name', data=P)
    f.close()