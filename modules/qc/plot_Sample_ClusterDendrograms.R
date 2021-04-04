qc.plotSample_ClusterDendrograms <- function(save=..saveImg_qc){
  dseq_qc <- qc.getDeseq()
  
  require(ggdendro)
  
  distances <- ..getCorrelationMatrix(getCounts(dseq_qc))
  labels <- dseq_qc$samplename
  
  hc <- hclust(distances, method = "complete")
  
  dend <- as.dendrogram(hc)
  dend_data <- dendro_data(dend, type = "rectangle")
  
  p <- ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE, labels = FALSE) + 
    theme(panel.background = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
    geom_text(data = dend_data$labels, aes(x, y, label = labels), hjust = 1, size = 3) +
    xlab("Hierarchical Clustering (Complete-linkage)") + ylab("Distance") +
    ylim(-max(distances)/4,max(distances)) +    # shows warning "Scale for 'y' is already present..."
    ggtitle(label = "Dendrogram of Sample Clusters")
  if(max(distances) > 0.1){
    p <- p + geom_hline(yintercept = 0.1, linetype="dashed", color="darkgray")
  }
  print(p)
  if(save) savePlot(p, 6, 6, filedir=qc.getResultDir(), filename="Cluster Dendrogram")
  
  invisible(p)
}
