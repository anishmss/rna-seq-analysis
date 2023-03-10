filter.genes <- function(dseq, keep_geneIds){
  dseq_filtered <- dseq[which(row.names(assay(dseq)) %in% unique(keep_geneIds)), ]
  return(invisible(dseq_filtered))
}

filter.samples <- function(dseq, removeSamples){
  dseq <- dseq[, !(dseq[['samplename']] %in% removeSamples)]
  dseq[['samplename']] <- droplevels(dseq[['samplename']])
  return(invisible(dseq))
}

filter.group <- function(dseq, removeLevel, fromFactor){
  dseq <- dseq[,dseq[[fromFactor]] != removeLevel]
  for(column in colnames(colData(dseq))){
    if(class(dseq[[column]]) == 'factor'){
      dseq[[column]] <- droplevels(dseq[[column]])
    }
  }
  return(invisible(dseq))
}

filter.printSummary <- function(dseq){
  writeLines(paste("Remaining gene count:", nrow(dseq)))
  writeLines(paste("Remaining sample count:", ncol(dseq)))
  
  writeLines("\nSample group counts:")
  cnttable <- data.frame(colData(dseq)) %>% 
    select("location", "condition") %>% 
    table() %>% 
    addmargins()
  print(cnttable)
  
  return(invisible())
}
