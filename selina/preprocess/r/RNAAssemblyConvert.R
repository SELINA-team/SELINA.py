RNAAssemblyConvert <- function(countMat, from = "GRCh37", to = "GRCh38",  dataPath) {
  countMat = as(as.matrix(countMat), "dgCMatrix")
  genes_original = rownames(countMat)
  load(paste0(dataPath, "/GRCh37.GRCh38.RData"))
  ensembl = GRCh37.GRCh38
  if (from == "GRCh37" & to == "GRCh38") {
    genes_convert = ensembl[match(genes_original, ensembl$Gene.name.GRCh37), "Gene.name.GRCh38"]
  }
  if (from == "GRCh38" & to == "GRCh37") {
    genes_convert = ensembl[match(genes_original, ensembl$Gene.name.GRCh38), "Gene.name.GRCh37"]
  }
  count_index = which((!duplicated(genes_convert)) & (!is.na(genes_convert)))
  countMat = countMat[count_index,]
  rownames(countMat) = genes_convert[count_index]
  return(countMat)
}
