# For plotting PCAs with a GUI tool
# BiocManager::install("pcaExplorer")

plotPCA_Sample_Interactive <- function(dseq_filtered){
  pcaExplorer::pcaExplorer(dseq_filtered, qc.getDeseq())
}
