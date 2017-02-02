#' Print model object
#' 
#' @param x A \code{\link{methimputeBinomialHMM}} object.
#' @param ... Ignored.
#' @return An invisible \code{NULL}.
#' @export
print.methimputeBinomialHMM <- function(x, ...) {
    
    print(x$data)
    message("Use the list operator $ to access all elements of this object.")
  
}
