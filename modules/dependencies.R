loadDependencies <- function(){
  # Data manipulation
  library(magrittr)
  library(tibble)
  library(dplyr)
  library(tidyr)
  # 
  # # DESeq2
  library(DESeq2)
  library(tximport)
  library(BiocParallel)
  # 
  # # Graphs
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  library(cowplot)
  # 
  # # Documentation
  library(docstring)
  library(beepr)
  
  # Interal modules
  # import("modules/utils/utils.R")
  # import("modules/utils/utils_graph.R")
  # import("modules/dataprep/dataprep.R")
  # import("modules/filter/filter.R")
  # import("modules/qualitycontrol/qualitycontrol.R")
  # import("modules/degtest/degtest.R")  
  DIR_MODULES <- getPath_fromProjWd("modules")
  moduleFiles <- list.files(path = DIR_MODULES, recursive = TRUE)
  for(m in moduleFiles){source(paste0(DIR_MODULES,"/", m))}
}
