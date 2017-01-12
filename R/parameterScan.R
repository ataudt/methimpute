#' Perform a parameter scan
#' 
#' Perform a parameter scan for an arbitrary parameter.
#' 
#' @param f A function for which to perform the scan.
#' @param param A character with the parameter for which to perform the scan.
#' @param values A vector with parameter values for which to perform the scan.
#' @param ... Other parameters passed through to \code{f}.
#' @return A data.frame with loglikelihood values.
parameterScan <- function(f, param, values, ...) {
  
    arglist <- list(values[1], ...)
    names(arglist)[1] <- param
    result <- list()
    for (val in values) {
        arglist[1] <- val
        model <- do.call(f, arglist)
        result[[as.character(val)]] <- data.frame(val, loglik=model$convergenceInfo$logliks[length(model$convergenceInfo$logliks)])
    }
    result <- do.call(rbind, result)
    names(result)[1] <- param
    return(result)
    
}