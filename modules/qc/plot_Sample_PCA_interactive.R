# For plotting other PCs
# BiocManager::install("pcaExplorer")

plotPCA_Sample_Interactive <- function(dseq_qc, dseq_filtered){
  pcaExplorer(dseq_filtered, dseq_qc)
}
