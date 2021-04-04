qc.plotSample_ClusterHeatmap <- function(factors = qc.getDesignFactors(), isDisplayVal=TRUE, save=..saveImg_qc){
  dseq_qc <- qc.getDeseq()
  levels_factor1 <- levels(dseq_qc[[factors[[1]]]])
  
  # 1. Plot all samples in one heatmap
  ..plotSample_ClusterHeatmap(dseq_qc, save = save)
  
  # 2. Plot samples grouped by first factor
  for(lvl in levels_factor1){
    ..plotSample_ClusterHeatmap(dseq_qc, subsetByFactor = factors[[1]], levelName = lvl, isDisplayVal = TRUE, save = save)
  }
  
  # 3. Plot all samples grouped by Location & Condition
  if(length(factors) == 2){
    dseq_qc[['..samplegrp']] <- factor(paste(dseq_qc[[factors[[1]]]], dseq_qc[[factors[[2]]]], sep = "_"))
    for(samplegrp in levels(dseq_qc[['..samplegrp']])){
      ..plotSample_ClusterHeatmap(dseq_qc, subsetByFactor = '..samplegrp', levelName = samplegrp, isDisplayVal = TRUE, save = save)
    }       
  }

  invisible(dseq_qc)
}

..plotSample_ClusterHeatmap <- function(dseq_qc, subsetByFactor = NA, levelName = NA, isDisplayVal = FALSE, save=..saveImg_qc){
  require(pheatmap)
  
  plotSubtitle <- ""
  if(!is.na(subsetByFactor)){
    dseq_qc <- dseq_qc[,dseq_qc[[subsetByFactor]] == levelName]
    plotSubtitle <- paste0(" (Subset - ", levelName,")")
  }
  
  distances <- ..getCorrelationMatrix(getCounts(dseq_qc))
  labels <- as.character(dseq_qc$samplename)
  
  mat_dist <- as.matrix(distances)
  rownames(mat_dist) <- paste(dseq_qc$condition, dseq_qc$location, sep=".")
  colnames(mat_dist) <- NULL
  
  # Code to make correlation colors consistent for all heatmap
  # 1. Create color palette of length 100
  colorPal <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(115)[1:100]
  
  # 2. Create 100 breaks from 0.8 to 1 (breaks is one element longer than color vector - ?pheatmap)
  #    When breaks do not cover the range of values, 
  #    then any value larger than max(breaks) will have the largest color 
  #    and any value lower than min(breaks) will get the lowest color.
  breaks <- seq(from=0.8, to=1, length.out = 101)
  
  p <- pheatmap(# Data
                1-mat_dist, 
                drop_levels = TRUE,
                
                # Elements
                clustering_distance_rows = distances,
                clustering_distance_cols = distances,
                
                # Colors
                breaks = breaks, # Used for mapping values to colors
                color = colorPal,
                border_color=NA,
                
                # Sizes
                fontsize_row = 8,
                fontsize = 6,
                fontsize_number = 14,
                
                # Labels
                labels_row = labels,
                legend = !isDisplayVal,
                display_numbers = isDisplayVal,
                number_format = "%.3f",
                main = paste0("Sample Correlation", plotSubtitle))
  print(p)
  if(save){
    filename <- paste0("Cluster w Correlation Matrix", plotSubtitle)
    savePlot(p, 8, 6, filedir=qc.getResultDir(), filename=filename)
  } 
}
