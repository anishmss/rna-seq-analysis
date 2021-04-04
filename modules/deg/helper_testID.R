..getShortenedCoeffname <- function(dseq, coeffName){
  factors <- all.vars(dseq@design)
  name_factor1 <- factors[[1]]
  name_factor2 <- factors[[2]]
  
  # Case 1 : coeffName is "locationBAT.conditionE"
  #          return "BAT.E"
  if(grepl(name_factor1, coeffName) && grepl(name_factor2, coeffName)){
    shortCoeffName <- gsub(name_factor1, "", coeffName) %>% 
              gsub(name_factor2, "", .)
    return(shortCoeffName)
  }
  
  # Case 2 : coeffName is "location_BAT_vs_CAG" or "condition_E_vs_C"
  #          return "BAT" or "E"
  if(grepl(name_factor1, coeffName) || grepl(name_factor2, coeffName)){
    shortCoeffName <- strsplit(coeffName, "_")[[1]][[2]]
    return(shortCoeffName)
  }
  return(coeffName)
}

..getFullCoeffName <- function(dseq, shortCoeffName){
  factors <- all.vars(dseq@design)
  
  name_factor1 <- factors[[1]]
  name_factor2 <- factors[[2]]
  levels_1 <- levels(dseq@colData[[name_factor1]])
  levels_2 <- levels(dseq@colData[[name_factor2]])
  
  if(!grepl("\\.", shortCoeffName)){
    # Case 1 : Not an interaction effect
    # shortCoeffName      coeffName
    # "BAT"               "location_BAT_vs_CAG"
    # "E"                 "condition_E_vs_C"
    if(shortCoeffName %in% levels_1){
      coeffName <- paste0(c(name_factor1,shortCoeffName,"vs",levels_1[[1]]), collapse = "_")
    }else if(shortCoeffName %in% levels_2){
      coeffName <- paste0(c(name_factor2,shortCoeffName,"vs",levels_2[[1]]), collapse = "_")
    }
  }else{
    # Case 2 : Is an interaction effect
    # shortCoeffName        coeffName
    # "BAT.E"               "locationBAT.conditionE"
    interactionEff <- strsplit(shortCoeffName, "\\.")[[1]]
    coeffName <- paste0(name_factor1, interactionEff[[1]], ".", name_factor2, interactionEff[[2]])
  }
  return(coeffName)
}

..getLevelsFromContrast <- function(dseq, contrast){
  # Given a contrast with interaction terms, this function returns the relevant levels in factor 1 (See example cases below)
  
  coeffs <- ..extractCoeffs_fromContrast(contrast)
  interactionEff <- sapply(coeffs, function(x) strsplit(x, "\\.")[[1]][[1]]) # "BAT.E" is split to "BAT"
  
  if(length(interactionEff) == 1){
    # Case 1 : one interaction effect
    # input contrast        :   "BAT.E"  
    # output levelNames     :   c("CAG", "BAT")
    name_factor1 <- all.vars(dseq@design)[[1]]
    refLevel_factor1 <- levels(dseq@colData[[name_factor1]])[[1]]
    levelNames <- c(refLevel_factor1, interactionEff[[1]])
  }else{
    # Case 2 : two interaction effects
    # input contrast        :   "BAT.E - BIC.E" 
    # output levelNames     :   c("BAT", "BIC")
    levelNames <- c(interactionEff[[1]], interactionEff[[2]])
  }
  return(levelNames)
}

..getInteractionCoeffs <- function(dseq, contrast){
  # Returns the left-handside (lhs) and right-handside (rhs) coefficents of a contrast
  # Useful for graphing interaction plots
  # Example input : BAT.E
  # Output:
  #     $lhs
  #     [1] "condition_E_vs_C"
  # 
  #     $rhs
  #     [1] "condition_E_vs_C"       "locationBAT.conditionE"
  coeffs <- ..extractCoeffs_fromContrast(contrast) 
  
  if(length(coeffs) == 1){
    coeffs <- sapply(c(contrast, strsplit(contrast, "\\.")[[1]][[2]]), function(x){..getFullCoeffName(dseq, x)}, USE.NAMES = FALSE)
    coeffs_lhs_rhs <- list(lhs = coeffs[2], 
                           rhs = c(coeffs[2], coeffs[1]))
  }else{
    coeffs_lhs_rhs <- lapply(coeffs, function(x){
      sapply(c(strsplit(x, "\\.")[[1]][[2]], x), function(x){..getFullCoeffName(dseq, x)}, USE.NAMES = FALSE)
    })
    names(coeffs_lhs_rhs) <- c("lhs", "rhs")
  }
  return(coeffs_lhs_rhs)
}

..getSamplesRelatedToContrast <- function(dseq, contrast){
  # Returns the names of the samples that are involved in the contrast
  # Useful for plotting
  #
  # Example 1
  # input contrast : "BAT"
  # output samples :  BAT_C_*, CAG_C_*
  #
  # Example 2 
  # input contrast : "BAT.E"
  # output samples :  BAT_C_*, BAT_E_*, CAG_C_*, CAG_E_*
  #
  # Example 3
  # input contrast : "BAT.E - BIC.E"
  # output samples :  BAT_C_*, BAT_E_*, BIC_C_*, BIC_E_*
  
  factors <- all.vars(dseq@design)
  
  levels_1 <- levels(dseq@colData[[factors[[1]]]])
  levels_2 <- levels(dseq@colData[[factors[[2]]]])
  refLevel_factor1 <- levels_1[[1]]
  refLevel_factor2 <- levels_2[[1]]
  
  coeffs <- ..extractCoeffs_fromContrast(contrast)
  samples <- c()
  for(coeff in coeffs){
    if(coeff %in% levels_1){
      # Case 1 : "BAT"
      # Output : c("CAG_C","BAT_C")
      s <- c(refLevel_factor1, coeff)
      s <- paste0(s, "_", refLevel_factor2)
    }else if(coeff %in% levels_2){
      # Case 2 : "E"
      # Output : c("CAG_C","CAG_E")
      s <- c(refLevel_factor2, coeff)
      s <- paste0(refLevel_factor1, "_", s)
    }else{
      # Case 3 : "BAT.E"
      # Output : c("CAG_C","CAG_E", "BAT_C", "BAT_E")
      interactionEff <- strsplit(coeff, "\\.")[[1]]
      s <- rep(c(refLevel_factor1, interactionEff[[1]]), each = 2)
      s <- paste0(s, "_", rep(c(refLevel_factor2, interactionEff[[2]]), 2))
    }
    if(length(samples) == 0){
      samples <- s
    }else{
      # Example testID input : "BAT.E - BIC.E"
      # Samples from "BAT.E" : s1 = c("CAG_C","CAG_E", "BAT_C", "BAT_E")
      # Samples from "BIC.E" : s2 = c("CAG_C","CAG_E", "BIC_C", "BIC_E")
      # Final output : outersect(s1, s2) = c("BAT_C", "BAT_E", "BIC_C", "BIC_E")
      samples <- c(setdiff(samples, s), setdiff(s, samples))
    }
  }
  samplenames <- dseq@colData[[1]]
  samplenames <- samplenames[grepl(paste0(samples, collapse = "|"), samplenames)]
  return(samplenames)
}

..extractCoeffs_fromContrast <- function(contrast){
  coeffs <- gsub(" ", "", contrast) %>%
    strsplit("\\+|\\-", perl = TRUE)
  return(coeffs[[1]])
}

