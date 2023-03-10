getCounts <- function(dseq, isNormalized = TRUE){ 
  if(class(dseq) == "DESeqTransform" || !isNormalized){ 
    mat_counts <- assay(dseq)   # Use this for QC plots (i.e. returns variance stabilized transformed counts)
  }else{    
    mat_counts <- counts(dseq, normalized=TRUE) # Use this for anything else
  }
  return(mat_counts)
}

getCounts_asDf <- function(dseq, isNormalized = TRUE){ 
  df_counts <- getCounts(dseq, isNormalized) %>% data.frame() %>% 
    rownames_to_column(var="gene_id") %>% 
    as_tibble()
  colnames(df_counts)[-1] <- as.vector(dseq$samplename)
  return(df_counts)
}