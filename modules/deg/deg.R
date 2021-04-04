..category_order <- c("Uu", "uU", "Dd", "dD", "UD", "DU")
..dseq_deg <- NULL
..signif <- 0.1

deg.signifTest <- function(dseq, contrasts, signif = ..signif, verbose = TRUE){
  #' Performs multiple adjustment across contrasts and genes
  #'
  #' @param dseq_fitted A DESeq object (can already be fitted or not)
  #' @param contrasts character vector. specified in mappedCoefficients
  #' @param signif numerical. significance threshold
  #' @param verbose logical. If TRUE, then prints the summary for each test.
  #' 
  dseq <- DESeq(dseq)
  list_signifTest <- list()
  for(contrast in contrasts){
    list_signifTest[[contrast]] <- ..signifTest_single(dseq , contrast, signif)
  }
  
  list_signifTest <- ..adjustPvals_acrossContrasts(list_signifTest, signif, verbose = verbose)

  fitResults <- Class_FitResults(dseq = dseq, list_signifTest = list_signifTest)
  invisible(fitResults)
}

..signifTest_single <- function(dseq, contrast, signif = ..signif){
  .g(coeffName, testContrast) %=% getTestContrast(dseq, contrast)
  
  res <- results(dseq, alpha = signif, contrast=testContrast)
  resShrunk <- lfcShrink(dseq, res=res, type = "ashr", quiet = TRUE) # The default type="apeglm" does not allow use of contrast
  
  sigTestRes <- Class_SignifTest(contrast = contrast, res = resShrunk)
  return(sigTestRes)
}

getTestContrast <- function(dseq, contrast){
  # input       => output
  # "BAT-BIC"   => list(c("locationBAT", "locationBIC"), c(0,1,-1,0,0,0))
  # "BAT+BAT.E" => list(c("locationBAT", "locationBAT.conditionE"), c(0,1,0,0,1,0))
  # "BAT"       => list("locationBAT", list("locationBAT"))
  
  contrast_vector <-strsplit(gsub(" ", "", contrast), "\\+|\\-", perl = TRUE)[[1]]
  coeffName <- sapply(X= contrast_vector, FUN = function(x) {..getFullCoeffName(dseq, x)}, USE.NAMES = FALSE)
  
  testContrast <- rep(0, length(resultsNames(dseq)))   # initialize test contrast to c(0,0,0,0,0,0)
  indexCoeffs <- match(coeffName, resultsNames(dseq))
  testContrast[indexCoeffs] <- 1
  
  if(grepl("-", contrast, fixed = TRUE)){  # if contrast contains minus sign, set coefficient to -1
    testContrast[indexCoeffs[2]] <- -1
  }
  
  return(list(coeffName, testContrast))
}

..adjustPvals_acrossContrasts <- function(list_signifTest, signif, verbose = TRUE){
  # Step 1 - Initiate new adjusted pvalues
  cnt_genes <- list_signifTest[[1]]@res@nrows
  pvals <- rep(NA, length(list_signifTest) * cnt_genes)
  
  for(i in 1:length(list_signifTest)){
    df_res <- list_signifTest[[i]]@res
    index <- (cnt_genes * (i-1) + 1):(cnt_genes * i)
    
    # Step 2 - For each gene of each contrast, assign a unique index
    #        - This is just to safely reassign the new adjusted pvalue back to the 
    #          original gene and contrast
    df_res[['index']] <- index
    
    # Step 3 - For each index, get the original unadjusted p-value
    pvals[index] <- df_res$pvalue
    
    # Step 4 - Genes that were filtered out by DESeq (padj == NA) is not included
    #          in the adjustment (also not part of total count in Step 6)
    na_index <- df_res %>%  
      data.frame() %>%  
      filter(is.na(padj)) %>% 
      .[['index']]
    pvals[na_index] <- NA        # set these pvalues to NA
  }
  
  # Step 5 - Adjust p-values using Benjamini Hochberg method
  padj_multContrast <- p.adjust(pvals, method="BH", n = sum(!is.na(pvals)))
  
  for(i in 1:length(list_signifTest)){
    sigTestRes <- list_signifTest[[i]]
    index <- (cnt_genes * (i-1) + 1):(cnt_genes * i)
    
    # Step 6 - Assign adjusted p-values to the padj column of DESeq result 
    sigTestRes@res$padj <- padj_multContrast[index]
    
    # Step 7 - Extract DE genes (can be removed later) 
    sigTestRes@df_degenes <- extractDEgenes(sigTestRes@res, signif)
    list_signifTest[[i]] <- sigTestRes
    
    if(verbose){
      cnt_degenes <- sum(sigTestRes@res$padj < signif, na.rm = TRUE)
      writeLines(paste("\nResult of [",sigTestRes@contrast,"]"))
      writeLines(paste0("Signif. DE genes : ", cnt_degenes))
      summary(sigTestRes@res)
    }
  }
  invisible(list_signifTest)
}

extractDEgenes <- function(res, signif=..signif){
  r <- res %>% 
    data.frame() %>% 
    rownames_to_column('gene_id') %>% 
    filter(padj < signif)
  return(r)
}
