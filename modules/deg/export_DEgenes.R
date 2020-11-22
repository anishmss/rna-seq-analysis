deg.exportDEGenes <- function(fitted, df_annot = NA, outDirSuffix="", extractCols=c()){
  filedir <- deg.getResultDir(addSubdir = paste0("DE gene list", outDirSuffix))

  if(all(is.na(df_annot))){
    for(resData in fitted@list_signifTest){
      saveCSV(resData@df_degenes, filedir = filedir, filename = paste0("DEgenes_", resData@testId))
    }
  }else{
    df_annot <- df_annot %>% select(all_of(c("gene_id", extractCols)))
    for(resData in fitted@list_signifTest){
      df_degenes <- resData@df_degenes   
      df_merged <- merge(df_degenes, df_annot, by="gene_id", all.x=TRUE)
      saveCSV(df_merged, filedir = filedir, filename = paste0("DEgenes_", resData@testId))
    }
  }
  return(invisible(fitted))
}
