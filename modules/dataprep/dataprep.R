rsem2Deseq <- function(in_dir, df_sampleInfo, design, txIn=TRUE, txOut=FALSE, df_tx2gene=NA){
  library("DESeq2", quietly = TRUE, warn.conflicts = FALSE)
  library("tximport", quietly = TRUE, warn.conflicts = FALSE)
  
  samplenames <- df_sampleInfo[,1]
  
  if(txIn){
    file_counts <- file.path(in_dir, paste0(samplenames, ".isoforms.results"))    
  }else{
    file_counts <- file.path(in_dir, paste0(samplenames, ".genes.results"))
  }

  txi <- tximport(file_counts, type = "rsem", txIn = txIn, txOut = txOut, tx2gene = df_tx2gene)
  deseq <- DESeqDataSetFromTximport(txi, colData = df_sampleInfo, design = as.formula(design))
  return(deseq)
}
