#' Beta mirror distribution
#' 
#' Get values for the beta mirror distribution
#' 
#' @param x A vector of quantiles.
#' @param emissionParams A matrix with emission parameters.
#' @param weights A table giving the weights for the two components.
#' @importFrom stats dbeta
#' 
dbetaMirror <- function(x, emissionParams, weights) {
  
    e <- emissionParams
    y <- weights["unmethylated"] * stats::dbeta(x, e["unmethylated", "a"], e["unmethylated", "b"]) + weights["methylated"] * stats::dbeta(x, e["methylated", "a"], e["methylated", "b"])
    return(y)
}


#' Beta symmetric distribution
#' 
#' Get values for the beta symmetric distribution
#' 
#' @param x A vector of quantiles.
#' @param emissionParams A matrix with emission parameters.
#' @param weights A table giving the weights for the two components.
#' @importFrom stats dbeta
#' 
dbetaSymmetric <- function(x, emissionParams, weights) {
  
    e <- emissionParams
    y <- weights["heterozygous"] * stats::dbeta(x, e["heterozygous", "a"], e["heterozygous", "b"])
    return(y)
}


# ============================================================================
# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)
# ============================================================================
dnbinom.size <- function(mean, var) {
    return(mean^2 / (var - mean))
}

dnbinom.prob <- function(mean, var) {
    return(mean/var)
}

dnbinom.mean <- function(size, prob) {
    return(size/prob - size)
}

dnbinom.var <- function(size, prob) {
    return( (size - prob*size) / prob^2 )
}

