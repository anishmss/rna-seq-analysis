..getDf_plotGene_Distribution <- function(dseq){
  require(reshape2)
  
  counts <- getCounts_asDf(estimateSizeFactors(dseq), isNormalized = TRUE) %>% 
            column_to_rownames("gene_id") 
  
  df_plotData <- as.matrix(log(counts + 1, base = 2)) %>% 
                 melt(varnames = c("gene_id", "samplename"), value.name = "count") %>% 
                 merge(dseq@colData) %>% 
                 data.frame() %>% 
                 unite(col = "group", all.vars(dseq@design), sep = "_", remove = FALSE)
  
  # Columns produced for plotting
  #   samplename  : Factor w/ 25 levels "BAT_C_1","BAT_C_2",..: 
  #   gene_id       : Factor w/ 6898 levels "TRINITY_DN10022_c0_g1",..:
  #   count       : num    9.42 5.07 9.1 7.67 6 ...
  #   group       : chr    "BAT_C" "BAT_C" "BAT_C" "BAT_C" ...
  #   <factor_1>  : Factor w/ 3 levels "CAG","BAT","BIC": 
  #   <factor_2>  : Factor w/ 2 levels "C","E": 
  
  return(df_plotData)
}

plotGene_Distribution <- function(dseq, saveImg = FALSE){
  df_plotData <- ..getDf_plotGene_Distribution(dseq)
  factors <- all.vars(dseq@design)
  
  # Graph 1 : Boxplot
  p <- ggplot(df_plotData, aes_string(x = "samplename", y = "count", fill = factors[[1]])) + 
    geom_boxplot() + 
    facet_wrap(~group,nrow=1, scale="free") +
    labs(x = "", y = "log2(normalized read counts + 1)", title = "Distribution of Read Counts") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Graph 2 : Density
  p2 <- ggplot(df_plotData, aes(x = count, colour = samplename, fill = samplename)) +
    geom_density(alpha = 0.2, size = 1.25) + 
    facet_wrap(~group, nrow=1) +
    labs(x = "log2(normalized read counts+1)", title ="Distribution of Read Counts") +
    theme(legend.position = "none") 
  
  print(p)
  print(p2)
  if(saveImg){
     savePlot(p, 12, 7, filedir=qc.getResultDir(), filename="Box Plot of Read Counts")
     savePlot(p2, 12, 7, filedir=qc.getResultDir(), filename="Density Plot of Read Counts")
  }
  
  # Returns a list of plots 
  invisible(
    list(p_boxplot = p, 
         p_density = p2)
  )
}