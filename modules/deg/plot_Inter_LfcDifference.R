deg.plotInter_LfcDifference <- function(fitted, tests, topX = 10, save=..saveImg_deg, returnPlot = FALSE){
  plots <- list()
  for(testId in tests[['testIDs']]){
    plots[[length(plots) + 1]] <- ..plotInter_LfcDifference(fitted, testId, topX, save)
  }
  
  if(returnPlot) return(invisible(plots))
  return(invisible(fitted))
}

..plotInter_LfcDifference <- function(fitted, testId, topX = 10, save=..saveImg_deg) {
  if(!"category" %in% colnames(fitted@list_signifTest[[testId]]@df_degenes)){
    fitted <- deg.addColumn_InterCategory(fitted, testId)
  }
  df_degenes <- fitted@list_signifTest[[testId]]@df_degenes
  
  df_interactionCtg <- df_degenes %>%
    mutate(difference = abs(lhs-rhs)) %>%
    group_by(category) %>%
    arrange(desc(difference)) %>%
    mutate(gene_id = factor(gene_id, levels=gene_id)) %>%    # levels tell ggplot not to sort gene_id's alphabetically, 
                                                             # but retain the original sort by descending `difference` 
    slice_head(n = topX) %>% 
    select(gene_id, lhs, rhs, category) %>%
    pivot_longer(cols = c("lhs", "rhs"), names_to="side", values_to="lfc", names_transform = list(side = as.factor))
  
  color_location <- brewer.pal(name="Dark2", n=3)
  color_location <- c("BAT" = color_location[1], "BIC" = color_location[2], "CAG" = color_location[3])
  locationsOfInterest <- unlist(..getInteractionGroups(testId))
  
  p <- ggplot(df_interactionCtg, aes(x = gene_id, y = lfc)) +
    geom_line(linetype="dashed", color="gray") +
    geom_abline(intercept = 0, slope = 0, color="gray20") +
    geom_point(aes(col=side)) +
    facet_grid(cols = vars(category), scales = "free_x", space = "free_x") +
    labs(title="Top 10 Log2 Fold Change difference per Category",
         subtitle=paste("Test:", testId),
         y="log2FoldChange",
         x="Genes") +
    scale_color_manual(name = "Location", labels = locationsOfInterest, values=as.vector(color_location[locationsOfInterest])) + 
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_blank()
    )
  
  print(p)
  if(save){
    savePlot(p, 9, 6, filedir=deg.getResultDir(addSubdir=paste0("Visualization/test - ", testId)), 
             filename="Top 10 genes per Interaction Category")
  }
  return(invisible(p))
}