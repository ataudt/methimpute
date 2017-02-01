#' methimpute objects
#'
#' @description
#' \pkg{\link{methimpute}} defines several objects.
#' \itemize{
#' \item \code{\link{methimputeData}}: Returned by \code{\link{importBSSeeker}}, \code{\link{importBismarck}} and \code{\link{inflateMethylome}}.
#' \item \code{\link{methimputeBinomialHMM}}: Returned by \code{\link{callMethylation}}.
#' }
#'
#' @name methimpute-objects
NULL


#' methimputeData
#'
#' A \code{\link[GenomicRanges]{GRanges}} object containing cytosine coordinates with meta-data columns 'context' and 'counts'.
#' @name methimputeData
#' @seealso \code{\link{methimpute-objects}}
NULL


#' methimputeBinomialHMM
#' 
#' The \code{methimputeBinomialHMM} is a list() which contains various entries (see Value section)
#' 
#' @return
#' A list() with the following entries:
#' \item{convergenceInfo}{A list() with information about the convergence of the model fitting procedure.}
#' \item{params}{A list() with fitted and non-fitted model parameters.}
#' \item{params.initial}{A list() with initial values for the model parameters.}
#' \item{data}{A \code{\link[GenomicRanges]{GRanges}} with cytosine positions and methylation status calls.}
#' \item{segments}{The \code{data} entry where coordinates of consecutive cytosines with the same methylation status have been merged.}
#' 
#' @name methimputeBinomialHMM
#' @seealso \code{\link{methimpute-objects}}
NULL

