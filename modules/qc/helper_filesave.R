..saveImg_qc <- FALSE

qc.setSaveImg <- function(saveImg){
  ..saveImg_qc <<- saveImg
}

qc.getResultDir <- function(addSubdir = ""){
  dir <- paste0("results/QC/", appendTimestamp("QC"))
  if(addSubdir != "")
    dir <- paste0(dir, "/", addSubdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  return(dir)
}
