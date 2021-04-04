deg.plot_BarGraph <- function(fitted, tests, save = ..saveImg_deg, returnPlot = FALSE){
  counts <- c()
  for(testId in tests[['contrasts']]){
    direction <- fitted@list_signifTest[[testId]]@df_degenes %>% 
      mutate(direction = as.factor(ifelse(log2FoldChange > 0, "up", "down"))) %>% 
      pull(direction)
    counts <- c(counts, sum(direction=="up"), sum(direction=="down"))
  }
  
  df_counts <- data.frame(labels = factor(rep(tests[['plotLabels']], each=2), levels=rev(tests[['plotLabels']])), 
                          direction = rep(c("up-regulated", "down-regulated"), length(tests[['contrasts']])),
                          counts = counts, 
                          stringsAsFactors = TRUE)
  
  p <- ggplot(data=df_counts, aes(x=counts, y=labels, fill=direction)) +
    geom_bar(stat="identity", position=position_dodge2(preserve = "single")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    xlab("DE gene counts") +
    ylab("Main Effects") +
    scale_fill_manual(values=c("#FF6962","#09BCA8")) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  print(p)
  
  if(save){
    savePlot(p, 8, 6, filedir=deg.getResultDir(addSubdir="Visualization"), filename="Bar graph of Main Effects")
  }
  if(returnPlot) return(invisible(p))
  return(invisible(fitted))
}
