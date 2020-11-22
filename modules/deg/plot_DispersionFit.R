deg.plot_DispersionFit <- function(fitted){
  plotDispEsts(fitted@dseq)
  return(invisible(fitted))
}
