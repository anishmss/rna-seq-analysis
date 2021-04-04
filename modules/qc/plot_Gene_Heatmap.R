library("pheatmap")
qc.plotGene_Heatmap <- function(topX=50, save=..saveImg_qc){
  dseq_qc <- qc.getDeseq()
  
  # Heatmap of count matrix
  sortOrder <- order(rowMeans(assay(dseq_qc)),decreasing=TRUE)[1:topX]
  df <- data.frame(colData(dseq_qc)) %>%  select("condition","location")
  
  p <- pheatmap(  # Data
                  assay(dseq_qc)[sortOrder,],
                  annotation_col=df,
                  
                  # Aesthetics
                  cluster_rows=FALSE, 
                  cluster_cols=FALSE, 
                  show_rownames=FALSE,
                  main = paste0("Top ", topX, " highly expressing 'genes'"))
  print(p)
  if(save) savePlot(p, 10, 8, filedir=qc.getResultDir(), filename="Heatmap Read Counts")
  
  return(invisible(p))
}
