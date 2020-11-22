getGeneId_toLabelMapping <- function(df_filter, columnAsLabel){
  map <- df_filter %>% pull(columnAsLabel)
  names(map) <- df_filter %>%  pull("gene_id")
  return(map)
}

getTop_Geneids <- function(df_degenes, topX=20, orderby="padj"){     # orderby : ["padj", "lfc"]
  if(orderby == "padj"){
    df_degenes <- df_degenes %>% arrange(padj)  
  }else if(orderby == "lfc"){
    df_degenes <- df_degenes %>% arrange(desc(abs(log2FoldChange)))
  }
  
  df_degenes <- df_degenes %>%  
    pull(gene_id) %>% 
    head(topX)
  return(df_degenes)
}


