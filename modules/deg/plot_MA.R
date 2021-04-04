deg.plot_MA <- function(fitted, tests, df_isDesc=NA, save = ..saveImg_deg, returnPlot = FALSE){
  plots <- list()
  
  y_max <- max(sapply(fitted@list_signifTest, function(x) quantile(abs(x@res$log2FoldChange), probs=0.995)[[1]]))
  
  for(i in 1:length(tests[['contrasts']])){
    plots[[length(plots)+1]] <- ..deg.plotMA_wDesc(fitted, tests[['contrasts']][i], tests[['plotLabels']][i], df_isDesc=df_isDesc, y_max=y_max, save=save)      
  }
  if(returnPlot) return(invisible(plots))
  return(invisible(fitted))
}

..deg.plotMA_wDesc <- function(fitted, testId, mainEffLabel, df_isDesc=NA, y_max = NA, save=..saveImg_deg){
  res <- fitted@list_signifTest[[testId]]@res
  if(all(is.na(df_isDesc))){
    df_isDesc <- data.frame(genes=row.names(res), isDescendant=FALSE)    
  } 
  
  plotData <- geneplotter::plotMA(res, returnData=TRUE) %>% 
    mutate(gene_id = row.names(res)) %>% 
    merge(df_isDesc, by="gene_id") %>% 
    mutate(isDEandDesc = factor(ifelse(isDE & isDescendant, "DE_Desc", ifelse(isDE, "DE", "insig")),
                                levels = c("DE_Desc", "DE", "insig"))) %>% 
    filter(mean > 0)

  if(is.na(y_max)){
    y_max <- quantile(abs(res$log2FoldChange), probs=0.995)[[1]]
  }
  y_min <- -y_max
  
  plotData_bounded <- plotData %>% mutate(outOfBounds = case_when(
    lfc > y_max ~  "top",
    lfc < y_min ~ "bottom",
    TRUE       ~ "no")) %>% 
    mutate(lfc = case_when(
      outOfBounds == "top" ~ y_max,
      outOfBounds == "bottom" ~ y_min,
      outOfBounds == "no" ~ lfc)) %>% 
    arrange(desc(isDEandDesc))
  
  
  p <- ggplot(plotData_bounded, aes(x=mean, y = lfc, col=isDEandDesc, shape=outOfBounds)) +
    geom_point(size=0.9) +
    geom_abline(intercept = 0, slope = 0, color = "gray30") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0.005,0.005), n.breaks = 10) +
    scale_color_manual(name = "", values = c("DE_Desc" = "red", "DE" = "black", "insig"="gray"), labels=c("DE_Desc" = "Stress Response DE Gene","DE" = "DE Gene","insig" = "Gene")) +
    scale_shape_manual(values = c(6, 16, 2), guide=FALSE) + 
    labs(title = paste0("MA plot of Treatment Effect for ", mainEffLabel), 
         subtitle = paste("Test:", testId), 
         y = "log2FoldChange", 
         x = "mean of normalized counts")  + 
    theme_test()
  print(p)
  if(save){
    savePlot(p, 7, 5, filedir=deg.getResultDir(addSubdir=paste0("Visualization/test - ", testId)),
                      filename="MA plot of Treatment Effect")
  }
  return(invisible(p))
}
