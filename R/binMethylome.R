#' Bin counts in windows
#' 
#' Bin counts from unmethylated and methylated cytosines in equidistant bins.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'counts.unmethylated' and 'counts.methylated'.
#' @param binsize The window size used for binning.
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#' 
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
#' @return A list() with \code{\link[GenomicRanges]{GRanges}} objects.
#' 
binMethylome <- function(data, binsize) {
  
    ptm <- startTimedMessage("Making fixed-width bins with ", binsize, "bp ...")
    chrom.lengths.floor <- floor(seqlengths(data) / binsize) * binsize
    tiles <- GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize)
    bins.template <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor, tilewidth=binsize), use.names=FALSE)
    seqlengths(bins.template) <- seqlengths(data)
    if (any(width(bins.template)!=binsize)) {
        stop("tileGenome failed")
    }
    stopTimedMessage(ptm)
    
    contexts <- c("total", levels(data$context))
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
        
        ptm <- startTimedMessage("Aggregating total methylation status and counts ...")
        data.context$status.methylated <- data.context$methylated == 1
        data.context$status.unmethylated <- !data.context$status.methylated
        ind <- findOverlaps(data.context, bins, select='first')
        dfa <- aggregate(as.data.frame(mcols(data.context)[,c('status.methylated', 'status.unmethylated', 'counts.unmethylated', 'counts.methylated')]), by=list(index=ind), FUN=sum, na.rm=TRUE)
        
        bins$status.methylated <- 0
        bins$status.methylated[dfa$index] <- dfa$status.methylated
        bins$status.unmethylated <- 0
        bins$status.unmethylated[dfa$index] <- dfa$status.unmethylated
        bins$status.ratio <- bins$status.methylated / (bins$status.methylated + bins$status.unmethylated)
        
        bins$counts.methylated <- 0
        bins$counts.methylated[dfa$index] <- dfa$counts.methylated
        bins$counts.unmethylated <- 0
        bins$counts.unmethylated[dfa$index] <- dfa$counts.unmethylated
        bins$counts.ratio <- bins$counts.methylated / (bins$counts.methylated + bins$counts.unmethylated)
        stopTimedMessage(ptm)
        
        bins.list[[context]] <- bins
    }
    
    return(bins.list)
}
