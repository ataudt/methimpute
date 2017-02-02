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

