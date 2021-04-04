Class_FitResults <- setClass("Class_FitResults", slots = list(
  dseq = "DESeqDataSet",
  list_signifTest = "list"
))

Class_SignifTest <- setClass("Class_SignifTest", slots = list(
  contrast = "character",
  df_degenes = "data.frame",
  res = "DESeqResults"
))
