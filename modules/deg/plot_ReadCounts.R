deg.plot_ReadCounts <- function(fitted, tests, labelFrom_df=NA, labelFrom_column=NA, topX = 20, orderby="padj", save=..saveImg_deg, returnPlot = FALSE){
  plots <- list()
  for(testId in tests[['testIDs']]){
    plots[[length(plots)+1]] <- ..plot_ReadCounts(fitted, testId, labelFrom_df=labelFrom_df, labelFrom_column=labelFrom_column, topX = topX, orderby=orderby, save=save)
  }
  if(returnPlot) return(invisible(plots))
  return(invisible(fitted))
}

..plot_ReadCounts <- function(fitted, testId, labelFrom_df=NA, labelFrom_column=NA, topX = 20, orderby="padj", save=..saveImg_deg){
  geneId_toLabel_map <- NA
  if(!all(is.na(labelFrom_df)) && !is.na(labelFrom_column)){
    geneId_toLabel_map <- getGeneId_toLabelMapping(labelFrom_df, labelFrom_column)
  }
  
  sigTestRes <- fitted@list_signifTest[[testId]]
  df_degenes <- sigTestRes@df_degenes
  
  data_plot <- getCounts_asDf(fitted@dseq, isNormalized = TRUE) %>% 
    filter(gene_id %in% getTop_Geneids(df_degenes, topX, orderby)) %>%
    select(c("gene_id", ..getSamplesOfInterest(colnames(.), testId))) %>%
    gather(colnames(.)[-1], key="samplename", value="counts") %>% 
    mutate(logcounts = log(counts + 1, base=2)) %>% 
    merge(., data.frame(colData(fitted@dseq)), by="samplename", all.x = TRUE)
   
  p <- ggplot(data_plot) +
    geom_point(aes(x = gene_id, y = logcounts, color = condition, shape = condition, group = factor(condition)),
               size = 3, alpha = 0.7) +
    # scale_y_log10() +
    xlab("Genes") +
    ylab("log2(counts + 1)") +
    ggtitle(label=paste0("Top ", topX, " Significant DE Genes (by ",orderby,")"), subtitle = paste("Test :", sigTestRes@testId)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  if(!all(is.na(geneId_toLabel_map))){
    p <- p + scale_x_discrete(labels=geneId_toLabel_map[data_plot$gene_id])
  }
  
  print(p) 
  if(save){
    savePlot(p, 10, 7, filedir=deg.getResultDir(addSubdir=paste0("Visualization/test - ", testId)), filename=paste0("Top DE Genes Read Count (top ",orderby,")"))
  }
  return(invisible(p))
}
