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

appendTimestamp <- function(filepath){
  if(!exists("..resultsTimestamp")) setTimestampNow()
    
  return( paste(filepath, ..resultsTimestamp) )
}

setTimestampNow <- function(){
  ..resultsTimestamp <<- format(Sys.time(), format = "(%Y-%m-%d_%H%M%S)")
  
  writeLines(paste0("Time set to ", ..resultsTimestamp))
}
