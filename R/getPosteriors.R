#' Get original posteriors
#' 
#' Transform the 'posteriorMeth', 'posteriorMax', and 'status' columns into original posteriors from the HMM.
#' 
#' @param data The \code{$data} entry from a \code{\link{methimputeBinomialHMM}} object.
#' @return A matrix with posteriors.
#' 
getPosteriors <- function(data) {
    
    p <- array(NA, dim = c(length(data), length(levels(data$status))), dimnames = list(NULL, status = levels(data$status)))
    
    mask <- data$status == 'Methylated'
    p[mask,'Methylated'] <- data$posteriorMax[mask]
    p[mask,'Intermediate'] <- (data$posteriorMeth[mask] - p[mask,'Methylated']) * 2
    p[mask,'Unmethylated'] <- 1 - p[mask,'Methylated'] - p[mask,'Intermediate']
    
    mask <- data$status == 'Intermediate'
    p[mask,'Intermediate'] <- data$posteriorMax[mask]
    p[mask,'Methylated'] <- (data$posteriorMeth[mask] - 0.5*p[mask,'Intermediate'])
    p[mask,'Unmethylated'] <- 1 - p[mask,'Intermediate'] - p[mask,'Methylated']
    
    mask <- data$status == 'Unmethylated'
    p[mask,'Unmethylated'] <- data$posteriorMax[mask]
    p[mask,'Intermediate'] <- (1 - p[mask,'Unmethylated'] - data$posteriorMeth[mask]) * 2
    p[mask,'Methylated'] <- 1 - p[mask,'Unmethylated'] - p[mask,'Intermediate']
    
    return(p)
}
