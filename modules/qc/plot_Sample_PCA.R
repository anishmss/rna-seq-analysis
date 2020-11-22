qc.plotSample_PCA <- function(dseq_qc, factors, showSamplenames=TRUE, save=..saveImg_qc, returnPlot=FALSE){  
  require(ggrepel)
  pcaData <- plotPCA(dseq_qc, intgroup=factors, ntop = length(dseq_qc), returnData=TRUE)

  if(length(factors) == 2){
    # Plot for Location AND Condition 
    p <- ..pcaData_toGGplot(pcaData) + 
      aes(color = dseq_qc[[factors[1]]], shape = dseq_qc[[factors[2]]]) +   
      labs(color=factors[1], shape=factors[2])
    
    if(length(levels(dseq_qc[[factors[2]]])) == 2){
      # Shapes for Condition : C = 1 (hollow circle); E = 0 (hollow square)
      p <- p + scale_shape_manual(values=c(1,0)) 
    }
  }else{
    # Plot for either Location OR Condition 
    p <- ..pcaData_toGGplot(pcaData) + 
      aes(color = dseq_qc[[factors]]) +
      labs(color=factors)
  }
  p <- p + ggtitle(label = "PCA of transformed counts")
  print(p)
  
  if(showSamplenames){
    pcaData <- pcaData %>% mutate(samplename = dseq_qc$samplename)
    p_labeled <-  p + geom_text_repel(data = pcaData, aes(label = samplename),
                                          box.padding   = 0.3,
                                          point.padding = 0.1,
                                          size = 3,
                                          segment.color = 'gray88')
    
    print(p_labeled)
  }
  
  if(save){
    savename <- paste0(factors, collapse = ".")
    savePlot(p, 6, 6, filedir=qc.getResultDir(), filename=paste0("PCA ", savename))

    if(showSamplenames)
      savePlot(p_labeled, 6, 6, filedir=qc.getResultDir(), filename=paste0("PCA (labeled) ", savename))
  }
  
  if(returnPlot){
    if(showSamplenames) return(p_labeled)
    return(p)
  }
  return(invisible(dseq_qc))
}

..pcaData_toGGplot <- function(pcaData){
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2)) +
    geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw()
  return(p)
}