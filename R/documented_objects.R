#' Toy data for arabidopsis (200.000bp of chr1)
#'
#' A \code{\link{methimputeData}} object for demonstration purposes in examples of package \pkg{\link{methimpute}}. The object contains the first 200.000 cytosines of chr1 from arabidopsis.
#'
#' @docType data
#' @name arabidopsis_toydata
#' @format A \code{\link{methimputeData}} object.
#' @examples
#'data(arabidopsis_toydata)
#'print(arabidopsis_toydata)
NULL

#' Gene coordinates for arabidopsis (chr1)
#'
#' A \code{\link[GenomicRanges]{GRanges}} object for demonstration purposes in examples of package \pkg{\link{methimpute}}. The object contains gene coordinates of chr1 from arabidopsis.
#'
#' @docType data
#' @name arabidopsis_genes
#' @format A \code{\link[GenomicRanges]{GRanges}} object.
#' @examples
#'data(arabidopsis_genes)
#'print(arabidopsis_genes)
NULL

#' Transposable element coordinates for arabidopsis (chr1)
#'
#' A \code{\link[GenomicRanges]{GRanges}} object for demonstration purposes in examples of package \pkg{\link{methimpute}}. The object contains transposable element coordinates of chr1 from arabidopsis.
#'
#' @docType data
#' @name arabidopsis_TEs
#' @format A \code{\link[GenomicRanges]{GRanges}} object.
#' @examples
#'data(arabidopsis_TEs)
#'print(arabidopsis_TEs)
NULL


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

