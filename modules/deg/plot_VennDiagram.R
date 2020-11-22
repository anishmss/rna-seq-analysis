deg.plot_VennDiagram <- function(fitted, tests, subtitle="", save = ..saveImg_deg, returnPlot = FALSE){
  list_degenes <- lapply(tests[['testIDs']], function(x) fitted@list_signifTest[[x]]@df_degenes[['gene_id']])
  names(list_degenes) <- tests[['labels']]

  require(purrr)
  require(RVenn)
  ..venn_object <- Venn(list_degenes)
  vd <- ggvenn(..venn_object,fill=c("gray5","gray50","gray100"))+
    coord_fixed(clip="off")+
    theme_void()+
    theme(legend.position = "none")
  
  if(save){ 
    savePlot(vd, 7, 7, filedir = deg.getResultDir(addSubdir = "Visualization"), 
             filename = paste("Venn Diagram", subtitle))
  }
  if(returnPlot) return(invisible(vd))
  return(invisible(fitted))
}
