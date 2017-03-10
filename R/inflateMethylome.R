#' Inflate an imported methylation extractor file
#' 
#' Inflate an imported methylation extractor file to contain all cytosine positions. This is useful to obtain a full methylome, including non-covered cytosines, because most methylation extractor programs only report covered cytosines.
#' 
#' @param methylome A \code{\link[GenomicRanges]{GRanges}} with methylation counts. 
#' @param methylome.full A \code{\link[GenomicRanges]{GRanges}} with positions for all cytosines or a file with such an object.
#' @return The \code{methylome.full} object with added metadata column 'counts'.
#' 
#' @export
#' @examples
#'## Get an example file in BSSeeker format
#'file <- system.file("extdata","arabidopsis_bsseeker.txt.gz", package="methimpute")
#'bsseeker.data <- importBSSeeker(file)
#'bsseeker.data
#'
#'## Inflate to full methylome (including non-covered sites)
#'data(arabidopsis_toydata)
#'full.methylome <- inflateMethylome(bsseeker.data, arabidopsis_toydata)
#'full.methylome
#'
inflateMethylome <- function(methylome, methylome.full) {
    
    if (is.character(methylome.full)) {
        ptm <- startTimedMessage("Loading full methylome from file ", methylome.full, " ...")
        temp.env <- new.env()
        methylome.full <- get(load(methylome.full, envir=temp.env), envir=temp.env) 
        stopTimedMessage(ptm)
    }
    if (is.character(methylome)) {
        ptm <- startTimedMessage("Loading methylome from file ", methylome, " ...")
        temp.env <- new.env()
        methylome <- get(load(methylome, envir=temp.env), envir=temp.env) 
        stopTimedMessage(ptm)
    }
    ptm <- startTimedMessage("Inflating methylome ...")
    counts <- array(0, dim=c(length(methylome.full), 2), dimnames=list(NULL, c("methylated", "total")))
    ind <- findOverlaps(methylome.full, methylome)
    counts[ind@from,] <- methylome$counts[ind@to,]
    methylome.full$counts <- counts
    stopTimedMessage(ptm)
    return(methylome.full)
    
}

