..getCorrelationMatrix <- function(x){
  return(as.dist(1 - cor(x, method = "pearson")))
}
