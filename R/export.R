#' Export a methylome
#' 
#' Export a methylome as a TSV file.
#' 
#' @param model A \code{\link{methimputeBinomialHMM}} object.
#' @param filename The name of the file to be exported.
#' @return \code{NULL}
#' @importFrom data.table fwrite
#' @export
#' @examples
#'## Get some toy data
#'file <- system.file("data","arabidopsis_toydata.RData", package="methimpute")
#'data <- get(load(file))
#'print(data)
#'model <- callMethylation(data, max.iter=10)
#'exportMethylome(model, filename = tempfile())
#'
exportMethylome <- function(model, filename) {

    ptm <- startTimedMessage("Writing to file ", filename, " ...")
    data <- model$data
    df <- methods::as(data, 'data.frame')
    df <- df[,c('seqnames', 'start', 'strand', 'context', 'counts.methylated', 'counts.total', 'posteriorMax', 'posteriorMeth', 'posteriorUnmeth', 'status','rc.meth.lvl')]
    data.table::fwrite(df,file = filename, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    stopTimedMessage(ptm)
    
}
