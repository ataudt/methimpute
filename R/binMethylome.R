#' Bin counts in windows
#' 
#' Bin counts from cytosines in equidistant bins.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'counts.unmethylated' and 'counts.methylated'.
#' @param binsize The window size used for binning.
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @importFrom stats aggregate
binCounts <- function(data, binsize) {
  
    ptm <- startTimedMessage("Making fixed-width bins with ", binsize, "bp ...")
    chrom.lengths.floor <- floor(seqlengths(data) / binsize) * binsize
    tiles <- GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize)
    bins <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize), use.names=FALSE)
    seqlengths(bins) <- seqlengths(data)
    if (any(width(bins)!=binsize)) {
        stop("tileGenome failed")
    }
    stopTimedMessage(ptm)
    
    ptm <- startTimedMessage("Aggregating counts ...")
    ind <- findOverlaps(data, bins, select='first')
    df <- stats::aggregate(as.data.frame(data$counts), by=list(ind), FUN=sum)
    bins$counts <- matrix(0, ncol=2, nrow=length(bins))
    colnames(bins$counts) <- colnames(data$counts)
    bins$counts[df$Group.1,'methylated'] <- df$methylated
    bins$counts[df$Group.1,'total'] <- df$total
    stopTimedMessage(ptm)
    
    return(bins)
}


#' Bin positions in windows
#' 
#' Bin cytosine positions in equidistant bins.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata column 'context'.
#' @param binsize The window size used for binning.
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#' 
binPositions <- function(data, binsize) {
  
    ptm <- startTimedMessage("Making fixed-width bins with ", binsize, "bp ...")
    chrom.lengths.floor <- floor(seqlengths(data) / binsize) * binsize
    tiles <- GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize)
    bins <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize), use.names=FALSE)
    seqlengths(bins) <- seqlengths(data)
    if (any(width(bins)!=binsize)) {
        stop("tileGenome failed")
    }
    stopTimedMessage(ptm)
    
    ptm <- startTimedMessage("Counting total overlaps ...")
    bins$total <- countOverlaps(bins, data)
    stopTimedMessage(ptm)
    
    ptm <- startTimedMessage("Counting context-specific overlaps ...")
    contexts <- levels(factor(data$context))
    bins$context <- matrix(0, ncol=length(contexts), nrow=length(bins), dimnames=list(NULL, context=contexts))
    for (context in contexts) {
        mask <- data$context == context
        bins$context[,context] <- countOverlaps(bins, data[mask])
    }
    stopTimedMessage(ptm)
    
    ptm <- startTimedMessage("Adding distance ...")
    bins$distance <- c(-1, start(bins)[-1] - end(bins)[-length(bins)] - 1)
    bins$distance[bins$distance < 0] <- Inf 
    stopTimedMessage(ptm)
    
    return(bins)
}


#' Bin cytosine positions in windows
#' 
#' Bin cytosine positions, methylation counts and state in equidistant windows.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'methylated', 'context', 'counts.unmethylated' and 'counts.methylated'.
#' @param binsize The window size used for binning.
#' @param contexts A character vector with contexts for which the binning will be done.
#' @param columns.average A character vector with names of columns in \code{data} that should be averaged in bins.
#' @return A list() with a \code{\link[GenomicRanges]{GRanges}} objects for each context.
#' 
#' @importFrom stats aggregate
binMethylome <- function(data, binsize, contexts='total', columns.average=c('posteriorMeth')) {
  
    ptm <- startTimedMessage("Making fixed-width bins with ", binsize, "bp ...")
    chrom.lengths.floor <- floor(seqlengths(data) / binsize) * binsize
    tiles <- GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize)
    bins.template <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize), use.names=FALSE)
    seqlengths(bins.template) <- seqlengths(data)
    if (any(width(bins.template)!=binsize)) {
        stop("tileGenome failed")
    }
    stopTimedMessage(ptm)
    
    bins.list <- list()
    for (context in contexts) {
        messageU("Working on context '", context, "'", underline=NULL, overline='-')
        bins <- bins.template
        if (context == 'total') {
            data.context <- data
        } else {
            mask <- data$context == context
            data.context <- data[mask]
        }
      
        ptm <- startTimedMessage("Counting cytosine overlaps ...")
        bins$cytosines <- countOverlaps(bins, data.context)
        stopTimedMessage(ptm)
        
        ptm <- startTimedMessage("Aggregating counts ...")
        ind <- findOverlaps(data, bins, select='first')
        df <- stats::aggregate(as.data.frame(data$counts), by=list(ind), FUN=sum)
        bins$counts <- matrix(0, ncol=2, nrow=length(bins))
        colnames(bins$counts) <- colnames(data$counts)
        bins$counts[df$Group.1,'methylated'] <- df$methylated
        bins$counts[df$Group.1,'total'] <- df$total
        stopTimedMessage(ptm)
        for (col in columns.average) {
            ptm <- startTimedMessage("Averaging column ", col, " ...")
            if (!is.null(mcols(data)[col])) {
                df <- stats::aggregate(mcols(data)[[col]], by=list(ind), FUN=mean)
                mcols(bins)[col] <- NA
                mcols(bins)[df$Group.1, col] <- df[,2]
            }
            stopTimedMessage(ptm)
        }
        
        bins.list[[context]] <- bins
    }
    
    return(bins.list)
}
