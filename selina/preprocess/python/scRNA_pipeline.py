import os
from selina.preprocess.python.scRNA_qc import scrna_qc


def preprocess_parser(subparsers):
    """
    Add argument parsers.
    """

    workflow = subparsers.add_parser(
        "preprocess",
        help="Run preprocess pipeline from scRNA-seq gene-cell count matrix.")

    group_input = workflow.add_argument_group("Input files arguments")
    group_input.add_argument("--format",
                             required=True,
                             dest="format",
                             default="",
                             choices=["h5", "mtx", "plain"],
                             help="Format of the count matrix file.")
    group_input.add_argument(
        "--matrix",
        dest="matrix",
        default="",
        help="Location of count matrix file. "
        "If the format is 'h5' or 'plain', users need to specify the name of the count matrix file."
        "If the format is 'mtx', the 'matrix' should be the name of .mtx formatted matrix file, such as 'matrix.mtx'."
    )
    group_input.add_argument(
        "--separator",
        dest="separator",
        default="tab",
        choices=["tab", "space", "comma"],
        help="The separating character (only for the format of 'plain')."
        "Values on each line of the plain matrix file will be separated by the character. DEFAULT: tab."
    )
    group_input.add_argument(
        "--feature",
        dest="feature",
        default="features.tsv",
        help="Location of feature file (required for the format of 'mtx'). "
        "Features correspond to row indices of count matrix. DEFAULT: features.tsv."
    )
    group_input.add_argument(
        "--gene-column",
        dest="gene_column",
        default=2,
        type=int,
        help=
        "If the format is 'mtx', please specify which column of the feature file to use for gene names. DEFAULT: 2."
    )
    group_input.add_argument(
        "--barcode",
        dest="barcode",
        default="barcodes.tsv",
        help="Location of barcode file (required for the format of 'mtx'). "
        "Cell barcodes correspond to column indices of count matrix. DEFAULT: barcodes.tsv. "
    )
    group_input.add_argument(
        "--gene-idtype",
        dest="gene_idtype",
        default="symbol",
        choices=["symbol", "ensembl"],
        help=
        "Type of gene name, 'symbol' for gene symbol and 'ensembl' for ensembl id. DEFAULT: symbol."
    )
    group_input.add_argument(
        "--assembly",
        dest="assembly",
        default="GRCh38",
        choices=["GRCh38", "GRCh37"],
        help="Assembly (GRCh38/hg38 and GRCh37/hg19). DEFAULT: GRCh38.")

    # Quality control cutoff
    group_cutoff = workflow.add_argument_group("Quality control arguments")
    group_cutoff.add_argument(
        "--count-cutoff",
        dest="count_cutoff",
        default=1000,
        type=int,
        help="Cutoff for the number of count in each cell. DEFAULT: 1000.")
    group_cutoff.add_argument(
        "--gene-cutoff",
        dest="gene_cutoff",
        default=500,
        type=int,
        help=
        "Cutoff for the number of genes included in each cell. DEFAULT: 500.")
    group_cutoff.add_argument(
        "--cell-cutoff",
        dest="cell_cutoff",
        default=10,
        type=int,
        help="Cutoff for the number of cells covered by each gene. DEFAULT: 10."
    )
    group_cutoff.add_argument(
        "--mito",
        dest="mito",
        action="store_true",
        help=
        "This flag should be used when you want to filter out cells with a high percentage of mitochondria genes."
    )
    group_cutoff.add_argument(
        "--mito-cutoff",
        dest="mito_cutoff",
        default=0.2,
        type=float,
        help=
        "Cutoff for the percentage of mitochondria genes in each cell. DEFAULT: 0.2."
    )

    group_process = workflow.add_argument_group("Process arguments")
    group_process.add_argument(
        "--variable-genes",
        dest="variable_genes",
        default=2000,
        type=int,
        help="Number of variable genes used in PCA. DEFAULT: 2000.")
    group_process.add_argument(
        "--npcs",
        dest="npcs",
        default=30,
        type=int,
        help="Number of dimensions after PCA. DEFAULT: 30.")
    group_process.add_argument("--cluster-res",
                               dest="cluster_res",
                               default=0.6,
                               type=float,
                               help="Clustering resolution. DEFAULT: 0.6.")

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument(
        "--directory",
        dest="directory",
        required=True,
        help=
        "Path to the directory where the result file shall be stored."
    )
    group_output.add_argument("--outprefix",
                              dest="outprefix",
                              default="query",
                              help="Prefix of output files. DEFAULT: query.")



# Generate Rscript
def GenerateRscript(count_file, gene_idtype, assembly, cell_cutoff, mito,
                    mito_cutoff, variable_genes, npcs, cluster_res, outprefix,
                    directory):

    rfile = os.path.join(directory, "%s.R" % (outprefix))
    outf = open(rfile, "w")
    refGenePath = os.path.split(os.path.split(
        os.path.abspath(__file__))[0])[0] + '/data'
    rsrcPath = os.path.split(os.path.split(
        os.path.abspath(__file__))[0])[0] + '/r'

    #========load package========
    script = '''# load package
        suppressMessages(library(Seurat))
        suppressMessages(library(ggplot2))
        suppressMessages(library(dplyr))
        suppressMessages(library(data.table))
        suppressMessages(library(presto))
        source("%s/RNARunSeurat.R")
        source('%s/RNAAssemblyConvert.R')
        source('%s/RNAEnsemblToSymbol.R')
        source('%s/FindMarkersMAESTRO.R')
    ''' % (rsrcPath, rsrcPath, rsrcPath, rsrcPath)
    outf.write(script)

    #========assembly conversion and gene id conversion========
    script = '''
        # read data
        expr = Read10X_h5("%s")
    ''' % (count_file)
    outf.write(script)

    if assembly == "GRCh37":
        if gene_idtype == "symbol":
            script = '''
                # assembly conversion
                expr = RNAAssemblyConvert(expr, from = "GRCh37", to = "GRCh38", dataPath = "%s")
            ''' % (refGenePath)
        elif gene_idtype == "ensembl":
            script = '''
                # gene id conversion
                expr = RNAEnsemblToSymbol(expr, assembly = "GRCh37", dataPath = "%s")
                expr = RNAAssemblyConvert(expr, from = "GRCh37", to = "GRCh38", dataPath = "%s")
            ''' % (refGenePath, refGenePath)
        outf.write(script)

    elif assembly == "GRCh38":
        if gene_idtype == "ensembl":
            script = '''
                # gene id conversion
                expr = RNAEnsemblToSymbol(expr, assembly = "GRCh38", dataPath = "%s")
            ''' % (refGenePath)
            outf.write(script)


#========analysis========
    if mito:
        mito = 'TRUE'
    else:
        mito = 'FALSE'

    script = '''
        # clustering
        RNA.res = RNARunSeurat(inputMat = expr, 
                            project = "%s", 
                            min.c = %d,
                            mito = as.logical("%s"),
                            mito.cutoff = %f,
                            variable.genes = %d, 
                            npcs = %d,
                            cluster.res = %f,
                            outdir = "%s",
                            outprefix = "%s")
    ''' % (outprefix, cell_cutoff, mito, mito_cutoff, variable_genes, npcs,
           cluster_res, directory, outprefix)
    outf.write(script)

    #========save seurat object========
    script = '''
        if (is.list(RNA.res)){
            # save object
            saveRDS(RNA.res[[1]], "%s_res.rds")
            # output differentially expressed genes
            write.table(RNA.res[[2]], file.path("%s", paste0("%s", "_cluster_DiffGenes.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
        } else {
            saveRDS(RNA.res, "%s_res.rds")
        }
    ''' % (os.path.join(directory, outprefix), directory, outprefix,
           os.path.join(directory, outprefix))
    outf.write(script)

    #========finish srcipt output========
    outf.close()
    return os.path.abspath(rfile)


def query_preprocess(directory, outprefix, fileformat, matrix, separator,
                     feature, gene_column, gene_idtype, barcode, count_cutoff,
                     gene_cutoff, cell_cutoff, mito, mito_cutoff,
                     variable_genes, npcs, cluster_res, assembly):

    try:
        os.makedirs(directory)
    except OSError:
        pass

    assembly = "GRCh38"
    scrna_qc(directory + "/Data", outprefix, fileformat, matrix, separator,
             feature, gene_column, barcode, count_cutoff, gene_cutoff,
             cell_cutoff, assembly)

    count_file = os.path.abspath(
        os.path.join(directory + "/Data",
                     outprefix + "_filtered_gene_count.h5"))

    rscript = GenerateRscript(count_file, gene_idtype, assembly, cell_cutoff,
                              mito, mito_cutoff, variable_genes, npcs,
                              cluster_res, outprefix, directory)

    cmd = "Rscript %s" % (rscript)
    os.system(cmd)
