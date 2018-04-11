#' Chromosome lengths for Arabidopsis
#'
#' A data.frame with chromosome lengths for Arabidopsis.
#'
#' @docType data
#' @name arabidopsis_chromosomes
#' @format A data.frame.
#' @examples
#'data(arabidopsis_chromosomes)
#'print(arabidopsis_chromosomes)
NULL

#' Toy data for Arabidopsis (200.000bp of chr1)
#'
#' A \code{\link{methimputeData}} object for demonstration purposes in examples of package \pkg{\link{methimpute}}. The object contains the first 200.000 cytosines of chr1 from Arabidopsis.
#'
#' @docType data
#' @name arabidopsis_toydata
#' @format A \code{\link{methimputeData}} object.
#' @examples
#'data(arabidopsis_toydata)
#'print(arabidopsis_toydata)
NULL

#' Gene coordinates for Arabidopsis (chr1)
#'
#' A \code{\link[GenomicRanges]{GRanges-class}} object for demonstration purposes in examples of package \pkg{\link{methimpute}}. The object contains gene coordinates of chr1 from Arabidopsis.
#'
#' @docType data
#' @name arabidopsis_genes
#' @format A \code{\link[GenomicRanges]{GRanges-class}} object.
#' @examples
#'data(arabidopsis_genes)
#'print(arabidopsis_genes)
NULL

#' Transposable element coordinates for Arabidopsis (chr1)
#'
#' A \code{\link[GenomicRanges]{GRanges-class}} object for demonstration purposes in examples of package \pkg{\link{methimpute}}. The object contains transposable element coordinates of chr1 from Arabidopsis.
#'
#' @docType data
#' @name arabidopsis_TEs
#' @format A \code{\link[GenomicRanges]{GRanges-class}} object.
#' @examples
#'data(arabidopsis_TEs)
#'print(arabidopsis_TEs)
NULL


#' methimpute objects
#'
#' @description
#' \pkg{\link{methimpute}} defines several objects.
#' \itemize{
#' \item \code{\link{methimputeData}}: Returned by \code{\link{importBSSeeker}}, \code{\link{importBismark}} and \code{\link{inflateMethylome}}.
#' \item \code{\link{methimputeBinomialHMM}}: Returned by \code{\link{callMethylation}}.
#' }
#'
#' @name methimpute-objects
NULL


#' methimputeData
#'
#' A \code{\link[GenomicRanges]{GRanges-class}} object containing cytosine coordinates with meta-data columns 'context' and 'counts'.
#' @name methimputeData
#' @seealso \code{\link{methimpute-objects}}
NULL


#' methimputeBinomialHMM
#' 
#' The \code{methimputeBinomialHMM} is a list() which contains various entries (see Value section). The main entry of this object is \code{$data}, which contains the methylation status calls and posterior values. See Details for a description of all columns.
#' 
#' The \code{$data} entry in this object contains the following columns:
#' \itemize{
#' \item{context }{The sequence context of the cytosine.}
#' \item{counts }{Counts for methylated and total number of reads at each position.}
#' \item{distance }{The distance in base-pairs from the previous to the current cytosine.}
#' \item{transitionContext }{Transition context in the form "previous-current".}
#' \item{posteriorMax }{Maximum posterior value of the methylation status call, can be interpreted as the confidence in the call.}
#' \item{posteriorMeth }{Posterior value of the "methylated" component.}
#' \item{posteriorUnmeth }{Posterior value of the "unmethylated" component.}
#' \item{status }{Methylation status.}
#' \item{rc.meth.lvl }{Recalibrated methylation level, calculated as \code{r$data$rc.meth.lvl = r$data$params$emissionParams$Unmethylated[data$context,] * r$data$posteriorUnmeth + r$params$emissionParams$Methylated[data$context,] * r$data$posteriorMeth}, where \code{r} is the \code{methimputeBinomialHMM} object.}
# #' \item{rc.counts }{Recalibrated counts for methylated and total number of reads at each position, calculated as \code{r$data$rc.counts[,2] <- sort(r$data$counts[,2])[rank(r$data$posteriorMax, ties.method = 'first')]; r$data$rc.counts[,1] <- round(r$data$rc.counts[,2] * r$data$rc.meth.lvl)}, where \code{r} is the \code{methimputeBinomialHMM} object.}
#' }
#' 
#' @return
#' A list() with the following entries:
#' \item{convergenceInfo}{A list() with information about the convergence of the model fitting procedure.}
#' \item{params}{A list() with fitted and non-fitted model parameters.}
#' \item{params.initial}{A list() with initial values for the model parameters.}
#' \item{data}{A \code{\link[GenomicRanges]{GRanges-class}} with cytosine positions and methylation status calls.}
#' \item{segments}{The \code{data} entry where coordinates of consecutive cytosines with the same methylation status have been merged.}
#' 
#' @name methimputeBinomialHMM
#' @seealso \code{\link{methimpute-objects}}
NULL

