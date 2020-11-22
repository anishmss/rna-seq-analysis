plotGene_ReadCount <- function(fitted, geneId, save=FALSE, returnPlot = FALSE){
  # Plot Of Normalized Counts For A Single Gene On Log Scale (joyce : but not log transformed)
  # Formula   : normalized counts plus a pseudocount of 0.5 are shown.
  # Reference : https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotCounts
  require(ggrepel)
  plotData <- plotCounts(fitted@dseq, gene=geneId, intgroup=c("condition", "location"), returnData = T)
  rownames(plotData) <- fitted@dseq$samplename
  
  p <- ggplot(plotData, aes(x = location, y = count, color = condition)) + 
    geom_point(position=position_jitter(w = 0.1,h = 0)) +
    ggtitle(label=paste0("Read counts of gene ", geneId)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- p + geom_text_repel(aes(label = rownames(plotData)))
  print(p)
  print(p2)
  if(save){
    savePlot(p, 10, 7, filedir=deg.getResultDir(addSubdir="Visualization"), filename=paste0("Read counts of gene ", geneId))
    savePlot(p2, 10, 7, filedir=deg.getResultDir(addSubdir="Visualization"), filename=paste0("Read counts of gene ", geneId, " (labeled)"))
  }
  if(returnPlot) return(invisible(p))
  return(invisible(fitted))
}
