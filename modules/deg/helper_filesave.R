..saveImg_deg <- FALSE

deg.setSaveImg <- function(saveImg){
  ..saveImg_deg <<- saveImg
}

deg.getResultDir <- function(addSubdir = ""){
  return(getResultDir(ifelse(addSubdir!= "",paste0("DE/", addSubdir), "DE")))
  # dir <- paste0("results/DE/", appendTimestamp("DE"))
  # if(addSubdir != "")
  #   dir <- paste0(dir, "/", addSubdir)
  # dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # return(dir)
}
