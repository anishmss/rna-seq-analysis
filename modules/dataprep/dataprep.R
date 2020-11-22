# tximport reference - http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#rsem

importRSEM <- function(dir_rsemCounts, fromProjRoot = TRUE, reload=FALSE){ # reload=TRUE will not read from .RDS file
  if(fromProjRoot)
    dir_rsemCounts <- getPath_fromProjWd(dir_rsemCounts)
  fileRDS <- paste0(dir_rsemCounts,".RDS")

  if(file.exists(fileRDS) && !reload){            # then load file, else, save it to fileRDS later
    writeLines(paste("Loading RSEM counts in ", fileRDS))
    dseq <- readRDS(fileRDS)
    return(dseq)
  }

  writeLines(paste("Loading RSEM counts in ", dir_rsemCounts))
  files <- list.files(dir_rsemCounts, pattern = ".genes.results")
  names(files) <- gsub(".genes.results","",files)
  txi <- tximport(paste0(dir_rsemCounts,"/",files), type="rsem", txIn=FALSE, txOut=FALSE)
  
  df_sampInfo <- data.frame(samplename = as.factor(names(files))) %>% 
    separate(samplename, into = c("location", "condition", "replicate"), sep = "_", remove=FALSE) %>% 
    mutate(location = factor(location, levels = FACTOR_LOC[["levels"]]),
           condition = factor(condition, levels = FACTOR_TRT[["levels"]]),
           samplegrp = factor(paste0(location, "_", condition), levels = ..getSubgroupLevels()))
  
  dseq <- DESeqDataSetFromTximport(txi,  colData = df_sampInfo, design = ~ location + condition + location:condition)
  dseq$condition <- relevel(dseq$condition, ref="C")
  dseq$location <- relevel(dseq$location, ref="CAG")
  if(!is.na(fileRDS) || reload){         # save dseq object to RDS for quicker loading
    saveRDS(dseq, fileRDS)
  }
  
  return(dseq)
}

addMetadata <- function(dseq, dir_metadata, fromProjRoot = TRUE){
  if(fromProjRoot)
    dir_metadata <- getPath_fromProjWd(dir_metadata)
  metadata <- read.csv(dir_metadata, stringsAsFactors = TRUE)
  
  df_sampInfo <- colData(dseq) %>% data.frame() %>% 
    inner_join(., metadata, by="samplename") %>% 
    mutate(weight_diff = weight_final - weight_init) %>% 
    mutate(carap_diff = carap_final - carap_init)
  
  for(col in colnames(df_sampInfo[,-c(1:3)])){
    dseq[[col]] <- df_sampInfo[,col]
  }
  writeLines("Metadata columns added : ")
  writeLines(colnames(colData(dseq))[-c(1:3)], sep = "    ")
  return(dseq)
}


addCategories <- function(dseq, categorize){
  toCategories <- function(df_sampInfo, colName, nbins = 5){
    column <- df_sampInfo[[colName]]
    if(all(is.na(column))){
      writeLines(paste(colName, "does not exist."))
      return(NA)
    }
    bin <- (max(column) - min(column)) / nbins
    ranges <- min(column) + (c(0:nbins) * bin)
    categorized <- cut(column, ranges, include.lowest = TRUE)
  }
  
  df_sampInfo <- colData(dseq) %>% data.frame()
  n <- length(categorize)
  
  for(i in 1:n){
    colName <- names(categorize[i])
    nbins <- categorize[[i]]
    out_colName <- paste0("ctg_",colName)
    
    dseq[[out_colName]] <- toCategories(df_sampInfo, colName, nbins)
    
    if(!all(is.na(dseq[[out_colName]]))){
      writeLines(paste0("Added categorized column: ",out_colName,"  (",nbins," bins)"))
    }
  }
  return(dseq)
}