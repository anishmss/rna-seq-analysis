deg.plotInter_Heatmap <- function(fitted, tests, transformMethod="vst", dseq_qc = NA, df_isDesc = NA, save=..saveImg_deg, returnPlot = FALSE){
  # Get variance-stabilized counts from either a prepared dseq_qc or by specifying the transformMethod (vst or rlog)
  if(!all(is.na(dseq_qc))){
    df_counts <- getCounts_asDf(dseq_qc)
  }else{
    df_counts <- getCounts_asDf(qcPrepare(fitted@dseq, method=transformMethod))
  }
  
  list_plots <- list()
  for(testId in tests[['testIDs']]){
    p <- ..plotInter_Heatmap(fitted, testId, df_counts, df_isDesc=df_isDesc, save = save)
    list_plots[[length(list_plots)+1]] <- p
  }
  
  if(returnPlot) return(invisible(list_plots))
  return(invisible(fitted))
}

..plotInter_Heatmap <- function(fitted, testId, df_counts, df_isDesc = NA, save=..saveImg_deg){
  if(!"category" %in% colnames(fitted@list_signifTest[[testId]]@df_degenes)){
    fitted <- deg.addColumn_InterCategory(fitted, testId)
  }
  df_degenes <- fitted@list_signifTest[[testId]]@df_degenes
  
  df_counts <- df_counts %>%
    filter(gene_id %in% df_degenes[["gene_id"]]) %>%
    select(c("gene_id", ..getSamplesOfInterest(dseq$samplename, testId))) %>% 
    merge(., select(df_degenes, gene_id, category), by="gene_id") %>% 
    arrange(category)
  
  require(ComplexHeatmap)
  # Transform to Z-score
  counts_Z <- df_counts %>% select(-one_of(c("gene_id", "category"))) %>%
    relocate(colnames(.)[grepl("CAG", colnames(.))])
  samplename <- colnames(counts_Z)
  counts_Z <- apply(counts_Z, 1, scale) %>% t()
  
  col_location <- brewer.pal(name="Dark2", n=3)
  col_location <- c("CAG" = col_location[3], "BAT" = col_location[1], "BIC" = col_location[2])
  top_annot = HeatmapAnnotation(
    location = factor(substr(samplename, 1,3), levels = unique(substr(samplename, 1,3))), 
    condition = as.factor(substr(samplename, 5,5)),
    col = list(location = col_location,
               condition = c("C" = "gray80", "E" = "gray60")),
    gp = gpar(col = "white"),
    show_annotation_name = c(location=FALSE, condition=FALSE),
    annotation_legend_param = list(location=list(nrow=1), condition=list(nrow=1)),
    annotation_label = list(location="site")
  )
  
  df_uniqueCtg <- df_counts %>% 
    select(category) %>% 
    mutate(index = row_number()) %>% 
    distinct(category, .keep_all = TRUE)
  left_annot <- rowAnnotation(
    category = anno_mark(at=df_uniqueCtg[,"index"], labels = df_uniqueCtg[,"category"],
                         side = "left")
  )
  marks <- character(nrow(df_counts))
  marks[df_uniqueCtg[["index"]]] <- as.character(df_uniqueCtg[["category"]])
  left_annot <- rowAnnotation(category = anno_text(marks, just = "right", location = 1, gp = gpar(fontsize=9)))
  library(circlize)
  col_fun <- colorRamp2(c(-max(abs(counts_Z)), 0, max(abs(counts_Z))), c("red", "black", "green"))
  
  if(!all(is.na(df_isDesc))){
    rightAnnot_data <- df_counts %>% select(gene_id) %>% 
      inner_join(., df_isDesc, by="gene_id") %>% 
      mutate(isDescendant = ifelse(isDescendant, 1, 0)) %>% 
      pull(isDescendant)
    right_annot <- rowAnnotation(stressResponse = anno_simple(rightAnnot_data, col = c("1"="black", '0'="white")), 
                                 show_annotation_name = c(stressResponse = FALSE))
  }else{
    right_annot <- NULL
  }
  
  h <- Heatmap(counts_Z, 
               col = col_fun,
               show_column_dend = FALSE,
               show_row_dend = FALSE,
               top_annotation = top_annot,
               # left_annotation = left_annot,
               right_annotation = right_annot,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               # row_split = df_counts$category,
               row_title = NULL,
               show_heatmap_legend = FALSE)
  draw(h)
  if(save){
    filedir <- deg.getResultDir(addSubdir=paste0("Visualization/test - ", testId))
    filename <- paste("Heatmap Count")
    finalPath <- paste0(filedir,"/",filename,".png")
    
    png(finalPath, width = 400, height = 600, units="px", res=100)
    draw(h)
    dev.off()
    writeLines(paste0("saving plot in ", finalPath))
  }
  return(invisible(h))
}

..getInteractionCtg <- function(lhs, rhs){
  if(any(is.na(c(lhs, rhs)))) return(NA)
  if(lhs > rhs){
    if(lhs > 0 && rhs >= 0) return("Uu")
    if(lhs <= 0 && rhs < 0) return("dD")
    if(lhs > 0 && rhs < 0) return("UD")
  }else{
    if(lhs >= 0 && rhs > 0) return("uU")
    if(lhs < 0 && rhs <= 0) return("Dd")
    if(lhs < 0 && rhs > 0) return("DU")
  }
}