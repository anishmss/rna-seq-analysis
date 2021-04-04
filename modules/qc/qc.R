..dseq_qc <- NULL
..design_qc <- NULL
..saveImg_qc <- FALSE

qc.normalize <- function(dseq, method="vst", blind=TRUE){  
  #' Variance Stabilizing Transformation
  #' 
  #' Transform counts for inducing homoskedasticity (similar variance)
  #' 
  if(method == "rlog"){
    # log norm w/ variance shrinkage
    ..dseq_qc <<- rlog(dseq, blind = blind)
  }else if(method == "vst"){
    ..dseq_qc <<- varianceStabilizingTransformation(dseq, blind = blind)
  }
  ..design_qc <<- dseq@design
  return(..dseq_qc)
}

qc.setSaveImg <- function(saveImg){
  ..saveImg_qc <<- saveImg
}

qc.getResultDir <- function(addSubdir = ""){
  return(getResultDir(ifelse(addSubdir!= "",paste0("QC/", addSubdir), "QC")))
  # dir <- paste0("results/QC/", appendTimestamp("QC"))
  # if(addSubdir != "")
  #   dir <- paste0(dir, "/", addSubdir)
  # dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # return(dir)
}


qc.setDeseq <- function(dseq){
  ..dseq_qc <<- dseq
}

qc.getDeseq <- function(){
  if(all(is.null(..dseq_qc))){
    stop("No DESeq2 object stored for QC. Call `qc.normalize()` first.")
  }
  return(..dseq_qc)
}

qc.getDesignFactors <- function(){
  return(all.vars(..design_qc))
}
