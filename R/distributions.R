#' Beta mirror distribution
#' 
#' Get values for the beta mirror distribution
#' 
#' @param emissionParams A matrix with emission parameters.
#' @param weights A table giving the weights for the two components.
#' 
dbetaMirror <- function(x, emissionParams, weights) {
  
    e <- emissionParams
    y <- weights["unmethylated"] * dbeta(x, e["unmethylated", "a"], e["unmethylated", "b"]) + weights["methylated"] * dbeta(x, e["methylated", "a"], e["methylated", "b"])
    return(y)
}


#' Beta symmetric distribution
#' 
#' Get values for the beta symmetric distribution
#' 
#' @param emissionParams A matrix with emission parameters.
#' @param weights A table giving the weights for the two components.
#' 
dbetaSymmetric <- function(x, emissionParams, weights) {
  
    e <- emissionParams
    y <- weights["heterozygous"] * dbeta(x, e["heterozygous", "a"], e["heterozygous", "b"])
    return(y)
}
