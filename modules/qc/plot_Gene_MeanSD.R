library(vsn)
qc.plotGene_MeanSD <- function(save=..saveImg_qc){
  dseq_qc <- qc.getDeseq()
  
  #' Plot transformed read counts for QC
  #' @description Assumption Check for Homoskedasticity
  
  plot <- meanSdPlot(assay(dseq_qc), ranks = F, plot = F)      # If ranks is False, data shows on the original scale
  plot$gg <- plot$gg +
    ggtitle(label="Mean SD of transformed read counts") +
    ylab("standard deviation")
  print(plot$gg)
  if(save)
    savePlot(plot$gg, 8, 7, filedir=qc.getResultDir(), filename="Mean SD")
  
  return(invisible(plot$gg))
}
