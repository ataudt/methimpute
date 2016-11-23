#' Make consensus segments
#' 
#' Given a list of \code{\link[GenomicRanges]{GRanges}}, construct a new GRanges as a \code{\link[GenomicRanges]{disjoin}} of the input segments.
#' 
makeConsensus <- function(models) {
    
  	## Load the files
  	models <- loadFromFiles(models, check.class='binnedMethylome')
  
  	## Get segments from list
  	ptm <- startTimedMessage("Getting segments ...")
  	segments.list <- GRangesList()
  	for (model in models) {
    		if (!is.null(model$segments)) {
      			segments.list[[as.character(model$ID)]] <- model$segments
    		}
  	}
  	stopTimedMessage(ptm)
  	
  	## Consensus GRanges
  	ptm <- startTimedMessage("Consensus ranges ...")
		consensus <- disjoin(unlist(segments.list, use.names=FALSE))
  	stopTimedMessage(ptm)
  	
  	## Consensus states
  	ptm <- startTimedMessage("Consensus states ...")
		constates <- matrix(NA, ncol=length(segments.list), nrow=length(consensus))
		for (i1 in 1:length(segments.list)) {
  			segments <- segments.list[[i1]]
  			splt <- split(segments, mcols(segments)$state)
  			mind <- as.matrix(findOverlaps(consensus, splt, select='first'))
  			# Assign Background state to NAs
  			mind[is.na(mind)] <- 1
  			constates[,i1] <- mind
		}
		consensus$states <- constates
  	stopTimedMessage(ptm)
		
# 		## Heterogeneity score
#     ptm <- startTimedMessage("Heterogeneity score ...")
#     S <- ncol(constates)
# 		physioState <- 0
#     tabs <- apply(constates - physioState, 1, function(x) { sort(table(x), decreasing=TRUE) }) # Heterogeneity score in reference to the physiological state
#     if (is.list(tabs) | is.vector(tabs)) {
#         consensus$Heterogeneity <- unlist(lapply(tabs, function(x) { sum(x * 0:(length(x)-1)) })) / S
#     } else if (is.matrix(tabs)) {
#         consensus$Heterogeneity <- colSums( tabs * 0:(nrow(tabs)-1) ) / S
#     }
#     stopTimedMessage(ptm)
  
  	## Variance
  	ptm <- startTimedMessage("Variance ...")
		consensus$variance <- apply(constates, 1, var)
    stopTimedMessage(ptm)
		
    return(consensus)
}