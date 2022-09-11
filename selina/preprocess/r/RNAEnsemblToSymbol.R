RNAEnsemblToSymbol <- function(countMat, assembly = "GRCh38", dataPath) {
  if (assembly == "GRCh38") {
    load(paste0(dataPath, '/GRCh38.ensembl.RData'))
    ensembl = GRCh38.ensembl
  } else if (assembly == "GRCh37") {
    load(paste0(dataPath, '/GRCh37.ensembl.RData'))
    ensembl = GRCh37.ensembl
  } 
  count_rowname = ensembl[match(rownames(countMat), ensembl$Gene.stable.ID), "Gene.name"]
  count_index = which(!duplicated(count_rowname) & !is.na(count_rowname))
  countMat = countMat[count_index,]
  rownames(countMat) = count_rowname[count_index]
  na_idx = which(is.na(count_rowname))
  if (length(na_idx) > 0) {
    warning(paste0("Omit ", length(na_idx), " genes of which symbol is not available !"))
  }
  return(countMat)
}