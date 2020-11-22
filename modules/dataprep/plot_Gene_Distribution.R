plotGene_Distribution <- function(dseq, normalized = TRUE, save = FALSE){
  require(reshape2)
  require(gridExtra)
  if(normalized){
    readCountType <- "normalized"    
    dseq <- estimateSizeFactors(dseq)
  }else{
    readCountType <- "raw"    
  }
  temp_data <- getCounts_asDf(dseq, isNormalized = normalized) %>% column_to_rownames("gene_id")
  temp_data <- as.matrix(log(temp_data + 1, base = 2))
  
  temp_df <- melt(temp_data)
  df <- data.frame(temp_df, Location = substr(temp_df$Var2, 1, 3), Condition = substr(temp_df$Var2, 5,5))
  df$samplegrp <- substr(df$Var2, 1,5)
  
  # Graph 1 : Boxplot
  p <-ggplot(df, aes(x = Var2, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
    ylab(paste0("log2(",readCountType," read counts + 1)")) +
    facet_wrap(~ samplegrp,nrow=1, scale="free") +
    ylim(min(df$value),max(df$value)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(label="Distribution of Read Counts")
  
  # Graph 2 : Density
  p2 <- ggplot(df, aes(x = value, colour = Var2, fill = Var2)) +
    geom_density(alpha = 0.2, size = 1.25) + 
    facet_wrap(~ samplegrp, nrow=1) +
    theme(legend.position = "none") + xlab(paste0("log2(",readCountType," read counts+1)")) +
    ggtitle(label="Distribution of Read Counts") 
  
  
  #gp <- grid.arrange(p, p2, ncol = 2,respect=TRUE) 
  #gp <- plot_grid(p,p2,ncol=2)
  
  print(p)
  print(p2)
  if(saveImg){
     savePlot(p, 10, 7, filedir=qc.getResultDir(), filename="Box Plot of Read Counts")
     savePlot(p2, 10, 7, filedir=qc.getResultDir(), filename="Density Plot of Read Counts")
  }
  invisible(dseq)
}