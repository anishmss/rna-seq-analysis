deg.plot_VennDiagram2 <- function(fitted, tests, subtitle="", save = ..saveImg_deg, returnPlot = FALSE,
                                 scale = 1, padding = 0.25, label_textSize = 7, count_textSize = 7, thickness = 1){
  padding <- padding * scale
  label_textSize  <- label_textSize * scale
  count_textSize <- count_textSize * scale
  thickness <- thickness * scale
  
  require(ggforce)
  require(limma)
  
  df_genes <- rowData(fitted@dseq) %>% data.frame() %>%  rownames_to_column("gene_id") %>% select("gene_id")
  for(testId in tests[['contrasts']]){
    df_genes[[testId]] <- df_genes$gene_id %in% fitted@list_signifTest[[testId]]@df_degenes[['gene_id']]
  }
  m_genes <- vennCounts(df_genes[,-1])
  class(m_genes)  <- "matrix"
  
  .g(df_counts, df_labels, df_shapes, circleFillColor) %=% ..getVennAxes(m_genes,tests, padding, length(tests[['contrasts']]))
  
  p <- ggplot(df_shapes) +
    geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .5, size = 1) +
    geom_circle(aes(x0 = x, y0 = y, r = 1.5), colour = 'black', lwd=thickness) +
    coord_fixed(clip = "off", ratio = 1) +
    theme_void() +
    theme(legend.position = 'None') +
    scale_fill_manual(values = circleFillColor) +
    labs(fill = NULL) +
    annotate("text", x = df_counts$x, y = df_counts$y, label = df_counts$Counts, size = count_textSize) +
    annotate("text", x = df_labels$x, y = df_labels$y, label = df_labels$label, size = label_textSize)
  print(p)
  
  if(save){  
    savePlot(p, width = 10, height = 7, filedir = deg.getResultDir(addSubdir="Visualization"), filename = paste0("Venn Diagram", subtitle))
  }
  if(returnPlot){
    return(invisible(p))
  }
  return(invisible(fitted))  
}


..getVennAxes <- function(m_genes, tests, padding, circles = 2){
  if(circles == 2){
    df_counts <- as.data.frame(as.matrix(m_genes)[-1,]) %>%
      mutate(x = c(1.33,  -1.33,   0),
             y = c(0.0,   0.0,   0)) %>% 
      select(x, y, Counts)
    
    df_labels <- data.frame(x=c(-1.33,1.33),
                            y=c(-1.5-padding, -1.5-padding),
                            label = tests[['plotLabels']])
    
    df_shapes <- data.frame(x = c(0.866, -0.866),
                            y = c(0, 0),
                            labels = tests[['contrasts']])
    circleFillColor <- c("gray50","gray100")
    
  }else if(circles == 3){
    df_counts <- as.data.frame(as.matrix(m_genes)[-1,]) %>%
      mutate(x = c(1.2,   -1.2,   0,    0,  0.8,  -0.8,   0),
             y = c(-0.6,  -0.6,  -1,  1.2,  0.5,   0.5,   0)) %>% 
      select(x, y, Counts)
    
    df_labels <- data.frame(x=c(0, -2.7-padding,2.7+padding),
                            y=c(2.8, -0.6-padding,-0.6-padding),
                            label = tests[['plotLabels']])
    
    df_shapes <- data.frame(x = c(0, 0.866, -0.866),
                            y = c(1, -0.5, -0.5),
                            labels = tests[['contrasts']])
    circleFillColor <- c("gray5","gray50","gray100")
  }else{
    stop('plot_VennDiagram2 can only plot 2-3 intersections')
  }
  return(list(df_counts=df_counts, df_labels=df_labels, df_shapes=df_shapes, circleFillColor=circleFillColor))
}
