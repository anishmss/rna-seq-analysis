library(kableExtra)

df_geneIDs <- group_by(getDfAnnot(), gene_id) %>% 
  summarise(swiss_prot = paste(swiss_prot, collapse = ", "), fly = paste(fly, collapse = ", "))

df_geneIsDesc <- read.csv(paste0(..PROJHOME, "/supplementary notebooks/outputs/geneID_isDescendant.csv"))      # Get this file from GO_stressResponse.Rmd


getDescendantGenes <- function(filtered = NA, geneids = NA){
  if(class(filtered) == "filteredDseq"){
    geneids <- row.names(counts(filtered@dseq))
  }
  df_descGenes <- df_geneIsDesc %>% 
    filter(gene_id %in% geneids) %>% 
    filter(isDescendant) 
  
  return(df_descGenes)
}

testOverrep <- function(total_genes, total_desc, cnt_degenes, cnt_desc){
  pval <- phyper(q = cnt_desc, m = total_desc, n = (total_genes - total_desc), k = cnt_degenes, lower.tail = FALSE) %>% round(4)
  # dhyper(x = 1:total_desc, m = total_desc, n = (total_genes - total_desc), k = cnt_degenes) %>% 
  #   plot() %>% 
  #   abline(v = cnt_desc)  
  return(pval)
}

trackDescGenes <- function(filtered, source="swiss_prot", fitted = NA, performTestIds = NA){  # source=c("swiss_prot","fly")
  if(class(fitted) != "modelFitDseq"){
    fitted <- modelfit(filtered)
  }
  
  df_descGenes <- getDescendantGenes(filtered = filtered)
  
  # Count how many are in outliers
  cooksCutoff <- qf(df1 = 6, df2 = length(fitted@dseq$samplename) - 6, p = .99)
  df_outliers <- mcols(fitted@dseq, use.names = TRUE) %>% 
    data.frame() %>% 
    subset(maxCooks > cooksCutoff) %>% 
    rownames_to_column("gene_id")
  
  cnt_totalGenes <- nrow(counts(fitted@dseq))
  df_outliersDesc <- df_descGenes %>%  filter(gene_id %in% df_outliers$gene_id)
  
  # writeLines(paste0("\nAfter filtering, the total number of genes in ",source," ID is : [ ",cnt_totalGenes, " ]"))
  # writeLines(paste0("After filtering, the number of descendant genes in ",source," ID is : [ ", nrow(df_descGenes), " ]"))
  # writeLines(paste0("After fitting, number of outliers are : [ ", nrow(df_outliers), " ]"))
  # writeLines(paste0("Number of descendants considered as an outlier : [ ", nrow(df_outliersDesc), " ]"))
  
  table_count <- c("Total number of genes" = cnt_totalGenes, 
                   "After filtering, the number of descendant genes is" = paste0(nrow(df_descGenes), " Desc"),
                   "After fitting, the number of outliers is" = paste0(nrow(df_outliersDesc)," Desc <br> out of ",nrow(df_outliers)))
  data.frame(
    "label" = names(table_count),
    "counts" = table_count, row.names = NULL) %>% 
    mutate( counts = cell_spec(counts, "html", color = ifelse(.[['counts']][2] == counts, "blue", 
                                                           ifelse(.[['counts']][3] == counts, "red", "default")))) %>% 
    kable(format = "html", escape = F,col.names = NULL) %>%
    kable_styling() %>% 
    writeLines()

  # Perform Hypothesis Testing
  if(!is.na(performTestIds)){
    fitted <- signifTest(fitted, list("BAT.E","BIC.E","BAT.E - BIC.E", "E", "E+BAT.E", "E+BIC.E"), verbose = FALSE)
  }
  list_DescCounts <- list()
  
  for(i in 1:length(fitted@list_plotResData)){
    resData <- fitted@list_plotResData[[i]]
    # Count how many are removed by independent filtering
    meanCutoff <- metadata(resData@res)$filterThreshold[[1]]
    
    df_lowMean <- data.frame(resData@res) %>% 
      rownames_to_column("gene_id") %>% 
      subset(!is.na(pvalue)) %>%
      subset(!gene_id %in% df_outliers$gene_id) %>%  # make sure outliers are already filtered out
      subset(baseMean < meanCutoff)
    
    df_lowMeanDesc <- df_descGenes %>%  filter(gene_id %in% df_lowMean$gene_id)
    
    # Count how many Descendants are left for testing
    cnt_descForTesting <- nrow(df_descGenes) - nrow(df_outliersDesc) - nrow(df_lowMeanDesc)
    
    # Count how many are DE genes
    list_degenes <- resData@df_degenes$genes
    resData@df_degenes[['isDescendant']] <- list_degenes %in% df_descGenes$gene_id
    cnt_DescDEgenes <- sum(resData@df_degenes[['isDescendant']])  
    
    # Save column isDescendant in df_degenes
    fitted@list_plotResData[[i]]@df_degenes <- resData@df_degenes
    
    overrepPval <- testOverrep(total_genes = cnt_totalGenes, total_desc = nrow(df_descGenes), 
                               cnt_degenes = length(list_degenes), cnt_desc = cnt_DescDEgenes)
    # writeLines(paste0("\nIn ", resData@testId))
    # writeLines(paste0("        Genes removed by independent filtering : ", nrow(df_lowMean), " (", nrow(df_lowMeanDesc), " Desc)"))
    # writeLines(paste0("        Final number of descendant genes remaining for testing : (", cnt_descForTesting, " Desc)"))
    # writeLines(paste0("        Number of DE genes : ",length(list_degenes), " (", cnt_DescDEgenes, " Desc)"))
    # writeLines(paste0("        Overrepresented Test p-value : ", overrepPval))
    
    
    list_DescCounts[[length(list_DescCounts) + 1]] <- list(testId = resData@testId, 
                                                           cnt_rmFilter = nrow(df_lowMean),
                                                           cnt_rmFilter_Desc = nrow(df_lowMeanDesc),
                                                           cnt_descForTesting = cnt_descForTesting,
                                                           cnt_degenes = length(list_degenes),
                                                           cnt_degenes_Desc = cnt_DescDEgenes,
                                                           overrepPval = overrepPval)
  }

  x <- data.frame(
    "Removed by Ind. Filtering" = apply(listOfList_toColumn(list_DescCounts, c("cnt_rmFilter_Desc", "cnt_rmFilter")), 2, function(x) paste0(x[[1]], " Desc out of ", x[[2]])),
    "Final descendants for testing" = paste0(listOfList_toColumn(list_DescCounts, "cnt_descForTesting"), " Desc"),
    "Signif. DE genes" = apply(listOfList_toColumn(list_DescCounts, c("cnt_degenes_Desc", "cnt_degenes")), 2, function(x) paste0(x[[1]], " Desc out of ", x[[2]])),
    "Overrepresented Test" = paste0("pval = ", listOfList_toColumn(list_DescCounts,"overrepPval"))
  , check.names = FALSE) %>% t()
  colnames(x) <- listOfList_toColumn(list_DescCounts, "testId")
  x  %>%  kable(format = "html", escape = F) %>%
    kable_styling() %>% 
    row_spec(3, bold = TRUE, color = "white", background = "#2B86FF") %>% 
    row_spec(4, italic = TRUE, color = "#AAA") %>% 
    writeLines()
  
  return(fitted)
}
