#' Bin counts in windows
#' 
#' Bin counts from unmethylated and methylated cytosines in equidistant bins.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'unmeth.counts' and 'meth.counts'.
#' @param binsize The window size used for binning.
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#' 
binCounts <- function(data, binsize) {
  
    ptm <- startTimedMessage("Making fixed-width bins ...")
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
    df <- aggregate(as.data.frame(mcols(data)[,c('unmeth.counts', 'meth.counts')]), by=list(ind), FUN=sum)
    bins$unmeth.counts <- 0
    bins$meth.counts <- 0
    bins$unmeth.counts[df$Group.1] <- df$unmeth.counts
    bins$meth.counts[df$Group.1] <- df$meth.counts
    stopTimedMessage(ptm)
    
    ptm <- startTimedMessage("Adding ratio and distance ...")
    bins$ratio <- bins$meth.counts / (bins$meth.counts + bins$unmeth.counts)
    bins <- bins[!is.na(bins$ratio)]
    bins$distance <- c(-1, start(bins)[-1] - end(bins)[-length(bins)] - 1)
    bins$distance[bins$distance < 0] <- Inf 
    stopTimedMessage(ptm)
    
    return(bins)
}