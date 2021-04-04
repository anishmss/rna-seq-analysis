savePlot <- function(plot, width, height, filedir, filename){
  dir.create(filedir, recursive=TRUE, showWarnings = FALSE)
  
  finalPath <- paste0(filedir,"/", filename,".png")
  print(paste0("saving plot in ", finalPath))
  ggsave(filename = finalPath, plot = plot, device = "png", width = width, height = height, units = c("in"))
}

saveCSV <- function(df, filedir, filename){
  dir.create(filedir, recursive=TRUE, showWarnings = FALSE)
  
  finalPath <- paste0(filedir,"/", filename,".csv")
  print(paste0("saving CSV in ", finalPath))
  write.csv(df, finalPath, row.names = FALSE)
}

saveParameters <- function(filename, ...){
  filename0 <- file.path(getResultDir(), filename)
  parameters <- lapply(substitute(list(...))[-1], deparse)
  parameters <- cbind(parameters, list(...))
  
  sink(filename0)
  print(parameters)
  sink()
  
  writeLines(paste0("saved in ", filename0))
}

getResultDir <- function(addSubdir = ""){
  dir <- file.path("results", getTimestamp())
  if(addSubdir != "") dir <- file.path(dir, addSubdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  return(dir)
}

getTimestamp <- function(){
  if(!exists("..resultsTimestamp")) setTimestampNow()
  
  return(..resultsTimestamp)
}

appendTimestamp <- function(filepath){
  return(paste0(filepath," (", getTimestamp(), ")") )
}

setTimestampNow <- function(){
  ..resultsTimestamp <<- format(Sys.time(), format = "%Y_%m_%d_%H%M")
  
  writeLines(paste0("Time set to ", ..resultsTimestamp))
}
