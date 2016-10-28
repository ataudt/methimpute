#' Bin counts in windows
#' 
#' Bin counts from unmethylated and methylated cytosines in equidistant bins.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'counts.unmethylated' and 'counts.methylated'.
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
    df <- aggregate(as.data.frame(mcols(data)[,c('counts.unmethylated', 'counts.methylated')]), by=list(ind), FUN=sum)
    bins$counts.unmethylated <- 0
    bins$counts.methylated <- 0
    bins$counts.unmethylated[df$Group.1] <- df$counts.unmethylated
    bins$counts.methylated[df$Group.1] <- df$counts.methylated
    stopTimedMessage(ptm)
    
    ptm <- startTimedMessage("Adding ratio and distance ...")
    bins$ratio <- bins$counts.methylated / (bins$counts.methylated + bins$counts.unmethylated)
    bins <- bins[!is.na(bins$ratio)]
    bins$distance <- c(-1, start(bins)[-1] - end(bins)[-length(bins)] - 1)
    bins$distance[bins$distance < 0] <- Inf 
    stopTimedMessage(ptm)
    
    return(bins)
}