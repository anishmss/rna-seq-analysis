deg.plotInter_CategoryTable <- function(fitted, tests, save=..saveImg_deg, returnPlot = FALSE){
  cnt_categories <- c()
  for(testId in tests[['testIDs']]){
    df_degenes <- fitted@list_signifTest[[testId]]@df_degenes
    if(!"category" %in% colnames(df_degenes)){
      fitted <- deg.addColumn_InterCategory(fitted, tests[['testIDs']])
    }
    
    cnt <- fitted@list_signifTest[[testId]]@df_degenes %>% 
      pull(category) %>% 
      table()
    cnt_categories <- c(cnt_categories, cnt)
  }
  
  df_cntCategory <- matrix(cnt_categories, byrow = TRUE, ncol = 6) %>% 
    data.frame() %>% 
    setNames(..category_order) %>% 
    mutate(label = tests[['labels']]) %>% 
    column_to_rownames("label")
  
  require(grid)
  require(gtable)
  require(gridExtra)
  p <- tableGrob(df_cntCategory, theme = ttheme_minimal()) %>% 
        gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)), t = 2, b = nrow(.), l = 1, r = ncol(.)) %>% 
        gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)), t = 1, l = 1, r = ncol(.))
  grid.newpage()
  grid.draw(p)
  
  if(save){
    filedir <- deg.getResultDir(addSubdir="Visualization")
    
    saveCSV(rownames_to_column(df_cntCategory, "category"), filedir = filedir, filename = "Counts per Interaction Category")
    savePlot(p, 4, 2, filedir = filedir, filename = "Counts per Interaction Category")
  }
  
  if(returnPlot) return(invisible(p))
  return(invisible(fitted))
} 
