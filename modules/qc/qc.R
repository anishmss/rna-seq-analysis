
qcPrepare <- function(dseq, method="vst", blind=TRUE){  
  #' Variance Stabilizing Transformation
  #' 
  #' Transform counts for inducing homoskedasticity (similar variance)
  #' 
  dseq_qc <- NA
  if(method == "rlog"){
    # log norm w/ variance shrinkage
    dseq_qc <- rlog(dseq, blind = blind)
  }else if(method == "vst"){
    dseq_qc <- varianceStabilizingTransformation(dseq, blind = blind)
  }
  
  return(dseq_qc)
}