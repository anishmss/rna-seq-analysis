# OPTION 1 : Plots several interaction plots in separate image files (shows coefficients table)
deg.plotInteractionPlot_Single <- function(fitted, tests, topX = 5, orderby = "lfc", labelFrom_df=NA, labelFrom_column=NA, df_isDesc = NA,  save = ..saveImg_deg, returnPlot = FALSE){
  for(testId in tests[['testIDs']]){
    sigTestRes <- fitted@list_signifTest[[testId]]
    df_degenes <- sigTestRes@df_degenes
    topGenes <- getTop_Geneids(df_degenes, topX = topX, orderby = orderby)
    subtitle <- paste0("(top ", orderby, ")")
    geneNameLabels <- ..getInterPlot_GeneNameLabel(topGenes, labelFrom_df, labelFrom_column, df_isDesc)
    
    plots <- list()
    for(i in 1:length(topGenes)){
      p <- ..plotInteractionPlot_Single(fitted, testId, geneName = topGenes[i], subtitle = subtitle, forMultiple=FALSE, geneNameLabel = geneNameLabels[i], save=save)
      plots[[length(plots)+1]] <- p
    }    
  }
  if(returnPlot) return(invisible(plots))
  return(invisible(fitted))
}


# OPTION 2 : Plots several interaction plots in 1 PDF file
deg.plotInteractionPlot_Multiple <- function(fitted, tests, topX=100, orderby="lfc", labelFrom_df=NA, labelFrom_column=NA, df_isDesc=NA){
  for(testId in tests[['testIDs']]){
    ..plotInteractionPlot_Multiple(fitted, testId, topX, orderby, labelFrom_df, labelFrom_column, df_isDesc)
  }
  return(invisible(fitted))
}


..plotInteractionPlot_Single <- function(fitted, testId=NA, geneName="", subtitle="", forMultiple = TRUE, geneNameLabel = NA, save=..saveImg_deg){
  require(ggpubr)
  require(gridExtra)
  sigTestRes <- fitted@list_signifTest[[testId]]
  df_degenes <- sigTestRes@df_degenes
  dseq <- fitted@dseq

  # Get gene's read counts
  df_counts <- getCounts_asDf(dseq)
  counts <- df_counts %>%
    filter(gene_id == geneName) %>% 
    select(colnames(.)[grepl(paste0(..getInteractionGroups(testId), collapse = "|"), colnames(.))]) %>% 
    gather(colnames(.)[-1], key = "samplename", value = "counts") %>% 
    merge(., data.frame(colData(dseq)), by="samplename") %>% 
    mutate(logcounts = log(counts + 1, base = 2))

  if(is.na(geneNameLabel))
    geneNameLabel <- geneName

  # Plot Graph
  p <- ggplot(counts, aes(x = condition, y = logcounts, color=location, shape=condition, group=location)) +
    geom_point(position=position_jitter(w = 0.1,h = 0), size=3) +
    ylab("log2(counts + 1)") +
    ggtitle(label = paste("Read counts of gene", geneNameLabel, subtitle), subtitle = paste0("Test: ", sigTestRes@testId)) +
    scale_shape_manual(values=c(1,0)) + 
    stat_summary(fun = mean, geom="line", size=1,  aes(color=location)) +
    theme_bw()
  
  
  if(forMultiple){
    # If for Multiple plots, remove legends and shrink texts
    p <- p + theme(legend.position = "None") + 
      labs(subtitle = geneNameLabel) + 
      theme(plot.subtitle = element_text(hjust = 0.5, size=7),
            axis.text.y.left = element_text(size=5),
            plot.title = element_blank(),
            axis.title = element_blank(),
            axis.text.x.bottom = element_blank(),
            axis.ticks.x = element_blank())
  }else{
    # If not forMultiple (i.e. plotting for one png file) show Coefficient Table
    
    coeffs <- data.frame(rowData(dseq))[geneName,] %>%  
      select("Intercept", unique(unlist(..getInteractionCoeffs(testId),use.names = FALSE))) %>% 
      t() %>% as.matrix() %>% data.frame() %>% 
      rownames_to_column("coefficient") %>% 
      dplyr::rename("estimate" = geneName)
    p2_theme <- ttheme_minimal(core = list(bg_params = list(fill = "gray95"),
                                           fg_params = list(hjust=0, x=0.1, fontsize=8)),
                               colhead = list(bg_params = list(fill = "gray95"),
                                              fg_params = list(fontsize=9, fontface="plain")))
    p2 <- tableGrob(coeffs, rows=NULL, theme = p2_theme)
    p <- ggarrange(p, p2, ncol = 2, nrow = 1, widths = c(3,1.5))
    print(p)
    
    if(save){
      savePlot(p, 8, 6, filedir=deg.getResultDir(addSubdir=paste0("Visualization/test - ", sigTestRes@testId)), filename=paste0("Interaction Plot ",subtitle," - ", geneName))
    }
  }
  return(invisible(p))
}

..plotInteractionPlot_Multiple <- function(fitted, testId, topX=100, orderby="lfc", labelFrom_df=NA, labelFrom_column=NA, df_isDesc=NA){
  sigTestRes <- fitted@list_signifTest[[testId]]
  df_degenes <- sigTestRes@df_degenes
  topGenes <- getTop_Geneids(df_degenes, topX = topX, orderby = orderby)
  subtitle <- paste0("(top ", orderby, ")")
  geneNameLabels <- ..getInterPlot_GeneNameLabel(topGenes, labelFrom_df, labelFrom_column, df_isDesc)
  
  pages <- ceiling(topX/20)
  for(page_i in 1:pages){
    topGenes_subset <- topGenes[((page_i-1)*20 + 1) : (page_i*20)] %>% na.omit()
    plots <- list()
    for(i in 1:length(topGenes_subset)){
      gene <- topGenes_subset[i]
      geneNameLabel <- geneNameLabels[((page_i-1)*20 + i)]
      
      p <- ..plotInteractionPlot_Single(fitted, testId, geneName = gene, subtitle = subtitle, forMultiple=TRUE, geneNameLabel = geneNameLabel, save=FALSE)
      plots[[length(plots) + 1]] <- p
    }    
    
    g <- ggarrange(plotlist = plots, common.legend = TRUE, legend = "right") %>% 
      arrangeGrob(., top=paste("Top DE genes of Test ",testId, "( by", orderby,")"))
    
    savePlot(g,  width=8.25, height = 9, 
             filedir = deg.getResultDir(addSubdir = paste0("Visualization/test - ", testId)),
             filename = paste0("Top Interaction Plots (by ",orderby,") (",page_i,")"))
  }
}

..getInterPlot_GeneNameLabel <- function(topGenes, labelFrom_df=NA, labelFrom_column=NA, df_isDesc=NA){
  if(!all(is.na(labelFrom_df)) && !is.na(labelFrom_column)){
    if(!all(is.na(df_isDesc))){
      df_isDesc$isDescendant <- ifelse(df_isDesc$isDescendant, "Y", "N")
      labelFrom_df <- merge(labelFrom_df, df_isDesc, by="gene_id", all.x=TRUE) %>% 
        mutate(!!labelFrom_column := paste0(.[[labelFrom_column]], " (", isDescendant, ")"))
    }
    geneId_toLabel_map <- getGeneId_toLabelMapping(labelFrom_df, labelFrom_column)
    geneNameLabels <- geneId_toLabel_map[topGenes]
  }else{
    geneNameLabels <- topGenes
  }
  return(geneNameLabels)
}