deg.addColumn_InterCategory <- function(fitted, testIDs){
  df_allGenes <- rowData(fitted@dseq) %>% 
    data.frame() %>%
    rownames_to_column("gene_id")
  
  for(testId in testIDs){
    .g(left_coeffs, right_coeffs) %=% ..getInteractionCoeffs(testId)
    df_degenes <- fitted@list_signifTest[[testId]]@df_degenes
    
    print(left_coeffs)
    print(right_coeffs)
    
    # This column is for filtering in Excel. Since it is not case insensitive, it sees Uu and uU as equal.
    categoryLabel <- paste0(1:6, "_", ..category_order)
    names(categoryLabel) <- ..category_order
    
    df_degenes <- df_allGenes %>% 
      filter(gene_id %in% df_degenes[["gene_id"]]) %>% 
      mutate(lhs = rowSums(.[left_coeffs]),
             rhs = rowSums(.[right_coeffs]),
             category = factor(mapply(..getInteractionCtg, lhs, rhs), levels = ..category_order),
             categoryLabel = categoryLabel[category]) %>% 
      select(gene_id, lhs, rhs, category, categoryLabel) %>% 
      merge(df_degenes, ., by="gene_id", all.x = TRUE)
    
    fitted@list_signifTest[[testId]]@df_degenes <- df_degenes
  }
  return(invisible(fitted))
}

  