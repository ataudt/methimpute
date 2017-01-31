#' Print model object
#' 
#' @param model A \code{\link{BinomialHMMcontext}} object.
#' @param ... Ignored.
#' @export
print.BinomialHMMcontext <- function(model, ...) {
    
    print(model$data)
    message("Use the list operator $ to access all elements of this object.")
  
}