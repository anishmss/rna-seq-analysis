FACTOR_LOC <- list(varName = "location", graphName = "site",
                   levels = c("CAG", "BAT", "BIC"))

FACTOR_TRT <- list(varName = "condition", graphName = "condition",
                   levels = c("C", "E"))

..category_order <- c("Uu", "uU", "Dd", "dD", "UD", "DU")

..getTestID_frCoeffName <- function(coeffName){
  varName_Loc <- FACTOR_LOC[['varName']]
  varName_Trt <- FACTOR_TRT[['varName']]
  
  # Case 1 : coeffName is "locationBAT.conditionE"
  #          return "BAT.E"
  if(grepl(varName_Loc, coeffName) && grepl(varName_Trt, coeffName)){
    testId <- gsub(varName_Loc, "", coeffName) %>% 
      gsub(varName_Trt, "", .)
    return(testId)
  }
  
  # Case 2 : coeffName is "location_BAT_vs_CAG" or "condition_E_vs_C"
  #          return "BAT" or "E"
  if(grepl(varName_Loc, coeffName) || grepl(varName_Trt, coeffName)){
    testId <- strsplit(coeffName, "_")[[1]][[2]]
    return(testId)
  }
  return(coeffName)
}


..getCoeffName_frTestID <- function(testID){
  # input         output
  # "BAT"       = "location_BAT_vs_CAG"
  # "BAT.E"     = "locationBAT.conditionE"
  # "E"         = "condition_E_vs_C")
  varName_Loc <- FACTOR_LOC[['varName']]
  varName_Trt <- FACTOR_TRT[['varName']]
  refLoc <- FACTOR_LOC[['levels']][[1]]
  refTrt <- FACTOR_TRT[['levels']][[1]]
  
  intEff <- strsplit(testID, "\\.")[[1]]
  if(length(intEff) == 1){
    if(testID %in% FACTOR_LOC[['levels']]){
      coeffName <- paste0(c(varName_Loc,testID,"vs",refLoc), collapse = "_")
    }else if(testID %in% FACTOR_TRT[['levels']]){
      coeffName <- paste0(c(varName_Trt,testID,"vs",refTrt), collapse = "_")
    }
  }else{
    coeffName <- paste0(varName_Loc, intEff[[1]], ".", varName_Trt, intEff[[2]])
  }
  return(coeffName)
}

..getInteractionGroups <- function(testID){
  testID_terms <- ..extractTestID_terms(testID)
  intEff <- sapply(testID_terms, function(x) strsplit(x, "\\.")[[1]][[1]])
  if(length(intEff) == 1){
    refLoc <- FACTOR_LOC[['levels']][[1]]
    samples <- c(refLoc, intEff[[1]])
  }else{
    samples <- c(intEff[[1]], intEff[[2]])
  }
  return(samples)
}

..getInteractionCoeffs <- function(testID){
  # input : BAT.E
  # [[1]]
  # [1] "condition_E_vs_C"
  # 
  # [[2]]
  # [1] "condition_E_vs_C"       "locationBAT.conditionE"
  testID_terms <- ..extractTestID_terms(testID)
  
  if(length(testID_terms) == 1){
    coeffs <- sapply(c(testID, strsplit(testID, "\\.")[[1]][[2]]), ..getCoeffName_frTestID, USE.NAMES = FALSE)
    coeffs_lhs_rhs <- list(lhs = coeffs[2], 
                           rhs = c(coeffs[2], coeffs[1]))
  }else{
    coeffs_lhs_rhs <- lapply(testID_terms, function(x){
      sapply(c(strsplit(x, "\\.")[[1]][[2]], x), ..getCoeffName_frTestID, USE.NAMES = FALSE)
    })
    names(coeffs_lhs_rhs) <- c("lhs", "rhs")
  }
  return(invisible(coeffs_lhs_rhs))
}

..getSamplesOfInterest <- function(samplenames, testID){
  searchPattern <- paste0(..getSampleGrpsOfInterest(testID), collapse = "|")
  samplenames <- samplenames[grepl(searchPattern, samplenames)]
  return(samplenames)
}

..getSampleGrpsOfInterest <- function(testID){
  refLoc <- FACTOR_LOC[['levels']][[1]]
  refTrt <- FACTOR_TRT[['levels']][[1]]
  testID_terms <- ..extractTestID_terms(testID)
  samples <- c()
  for(t in testID_terms){
    if(t %in% FACTOR_LOC[['levels']]){
      # Case 1 : "BAT"
      # Output : c("CAG_C","BAT_C")
      s <- c(refLoc, t)
      s <- paste0(s, "_", refTrt)
    }else if(t %in% FACTOR_TRT[['levels']]){
      # Case 2 : "E"
      # Output : c("CAG_C","CAG_E")
      s <- c(refTrt, t)
      s <- paste0(refLoc, "_", s)
    }else{
      # Case 3 : "BAT.E"
      # Output : c("CAG_C","CAG_E", "BAT_C", "BAT_E")
      intEff <- strsplit(t, "\\.")[[1]]
      s <- rep(c(refLoc, intEff[[1]]), each = 2)
      s <- paste0(s, "_", rep(c(refTrt, intEff[[2]]), 2))
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
  return(samples)
}

..getSubgroupLevels <- function(){
  levelsLoc <- FACTOR_LOC[['levels']]
  levelsTrt <- FACTOR_TRT[['levels']]
  l <- paste0(rep(levelsLoc, each=length(levelsTrt)), "_",
              rep(levelsTrt, length(levelsLoc)))
  return(l)
}

..extractTestID_terms <- function(testID){
  testID_terms <- gsub(" ", "", testID) %>%
    strsplit("\\+|\\-", perl = TRUE)
  return(testID_terms[[1]])
}

