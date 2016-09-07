#' Import a Bismarck methylation extractor file
#' 
#' Import a Bismarck methylation extractor file into a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param files The files to import.
#' @param chrom.lengths.df A data.frame obtained with \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for the seqlengths.
#' @param temp.store A folder where to save temporary files on disk. Set to \code{NULL} if no temporary files should be created.
#' @return A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'methylated' and 'context'.
#' 
importBismarck2 <- function(files, chrom.lengths.df, temp.store=tempfile("importBismarck")) {
  
  	classes <- c('NULL', 'character', 'character', 'numeric', 'character')
  	datas <- GRangesList()
  	
	  if (!is.null(temp.store)) {
        ## Create folder for temporary files
      	if (!file.exists(temp.store)) {
      	    dir.create(temp.store)
      	}
	  }
  	
  	### Do each file at a time to avoid memory overflow ###kk
  	data.dedups <- GRangesList()
  	for (i1 in 1:length(files)) {
  	    file <- files[i1]
	      savename <- file.path(temp.store, paste0(basename(file), '.RData'))
	      if (!file.exists(savename)) {
      	    ptm <- startTimedMessage("Reading file ", file, " ...")
        	  data.raw <- read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
          	data <- GRanges(seqnames=data.raw$V3, ranges=IRanges(start=data.raw$V4, end=data.raw$V4), strand="*", methylated=factor(data.raw$V2, levels=c("-", "+")), context=data.raw$V5)
          	rm(data.raw)
          	stopTimedMessage(ptm)
          	if (!is.null(temp.store)) {
          	    ptm <- startTimedMessage("Saving to file ", savename, " ...")
        	      save(data, file=savename)
              	stopTimedMessage(ptm)
          	}
	      } else {
      	    ptm <- startTimedMessage("Loading file ", savename, " ...")
      	    temp.env <- new.env()
      	    data <- get(load(savename, envir=temp.env), envir=temp.env)
          	stopTimedMessage(ptm)
	      }
	      
      	## Sorting
      	ptm <- startTimedMessage("Sorting ...")
      	data <- sort(data)
      	stopTimedMessage(ptm)
      	
  	  	## Get counts at each position
      	ptm <- startTimedMessage("Getting counts ...")
      	data$methylated <- as.logical(as.integer(data$methylated)-1)
      	data$mask.dup <- c(FALSE, as.logical(data@seqnames[-1] == data@seqnames[-length(data)])) & c(FALSE, as.logical(data@ranges@start[-1] == data@ranges@start[-length(data)]))
      	data$unmeth.sum <- rev(cumsum(rev(!data$methylated)))
      	data$meth.sum <- rev(cumsum(rev(data$methylated)))
      	data.dedup <- data[!data$mask.dup]
      	rm(data)
      	data.dedup$unmeth.counts <- -diff(c(data.dedup$unmeth.sum, 0))
      	data.dedup$meth.counts <- -diff(c(data.dedup$meth.sum, 0))
      	data.dedup$methylated <- NULL
      	data.dedup$mask.dup <- NULL
      	data.dedup$unmeth.sum <- NULL
      	data.dedup$meth.sum <- NULL
      	data.dedups[[length(data.dedups)+1]] <- data.dedup
      	stopTimedMessage(ptm)
      	
      	## Merging counts from data objects
      	if (length(data.dedups) == 2) {
          	ptm <- startTimedMessage("Merging counts ...")
          	data <- unlist(data.dedups, use.names=FALSE)
          	data <- sort(data)
          	rm(data.dedups)
          	data$mask.dup <- c(FALSE, as.logical(data@seqnames[-1] == data@seqnames[-length(data)])) & c(FALSE, as.logical(data@ranges@start[-1] == data@ranges@start[-length(data)]))
          	data$unmeth.sum <- rev(cumsum(rev(data$unmeth.counts)))
          	data$meth.sum <- rev(cumsum(rev(data$meth.counts)))
          	data.dedup <- data[!data$mask.dup]
          	rm(data)
          	data.dedup$unmeth.counts <- -diff(c(data.dedup$unmeth.sum, 0))
          	data.dedup$meth.counts <- -diff(c(data.dedup$meth.sum, 0))
          	data.dedup$mask.dup <- NULL
          	data.dedup$unmeth.sum <- NULL
          	data.dedup$meth.sum <- NULL
          	data.dedups <- GRangesList(data.dedup)
          	stopTimedMessage(ptm)
      	}
      	
  	}
  	data <- data.dedups[[1]]
  	rm(data.dedups, data.dedup)
  
  	## Assign seqlengths
  	chrom.lengths <- chrom.lengths.df$UCSC_seqlength
  	if (grepl('^chr',seqlevels(data)[1])) {
  		names(chrom.lengths) <- chrom.lengths.df$UCSC_seqlevel
  	} else {
  		names(chrom.lengths) <- chrom.lengths.df$NCBI_seqlevel
  	}
  	seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
  	
  	## Add ratio and distance
  	data$ratio <- data$meth.counts / (data$unmeth.counts + data$meth.counts)
    data$distance <- c(-1, start(data)[-1] - end(data)[-length(data)] - 1)
    data$distance[data$distance < 0] <- Inf 
  	return(data)
}


#' Import a Bismarck methylation extractor file
#' 
#' Import a Bismarck methylation extractor file into a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param files The files to import.
#' @param chrom.lengths.df A data.frame obtained with \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for the seqlengths.
#' @param temp.store A folder where to save temporary files on disk. Set to \code{TRUE} if memory is limiting.
#' @return A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'methylated' and 'context'.
#' 
importBismarck1 <- function(files, chrom.lengths.df, temp.store=NULL) {
  
  	classes <- c('NULL', 'character', 'character', 'numeric', 'character')
  	datas <- GRangesList()
  	
	  if (!is.null(temp.store)) {
        ## Create folder for temporary files
      	if (!file.exists(temp.store)) {
      	    dir.create(temp.store)
      	}
      	datafiles <- character()
        ## Import the files
      	for (i1 in 1:length(files)) {
      	    file <- files[i1]
    	      savename <- file.path(temp.store, paste0(basename(file), '.RData'))
    	      datafiles[[i1]] <- savename
    	      if (!file.exists(savename)) {
          	    ptm <- startTimedMessage("Reading file ", file, " ...")
            	  data.raw <- read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
              	data.gr <- GRanges(seqnames=data.raw$V3, ranges=IRanges(start=data.raw$V4, end=data.raw$V4), strand="*", methylated=factor(data.raw$V2, levels=c("-", "+")), context=data.raw$V5)
              	rm(data.raw)
              	stopTimedMessage(ptm)
          	    ptm <- startTimedMessage("Saving to file ", savename, " ...")
        	      save(data.gr, file=savename)
              	stopTimedMessage(ptm)
              	rm(data.gr)
    	      }
      	}
      	## Loading stored files
    	  for (i1 in 1:length(datafiles)) {
      	    ptm <- startTimedMessage("Loading file ", datafiles[[i1]], " ...")
      	    temp.env <- new.env()
      	    datas[[i1]] <- get(load(datafiles[[i1]], envir=temp.env), envir=temp.env)
          	stopTimedMessage(ptm)
    	  }
  	} else {
        ## Import the files
      	for (i1 in 1:length(files)) {
      	    file <- files[i1]
      	    ptm <- startTimedMessage("Reading file ", file, " ...")
        	  data.raw <- read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
          	data.gr <- GRanges(seqnames=data.raw$V3, ranges=IRanges(start=data.raw$V4, end=data.raw$V4), strand="*", methylated=factor(data.raw$V2, levels=c("-", "+")), context=data.raw$V5)
          	rm(data.raw)
          	stopTimedMessage(ptm)
    	      datas[[i1]] <- data.gr
      	}
      	rm(data.gr)
  	}
  	
  	## Sorting
  	ptm <- startTimedMessage("Sorting ...")
  	data <- unlist(datas, use.names=FALSE)
  	rm(datas)
  	data <- sort(data)
  	stopTimedMessage(ptm)
      	
  	## Assign seqlengths
  	chrom.lengths <- chrom.lengths.df$UCSC_seqlength
  	if (grepl('^chr',seqlevels(data)[1])) {
  		names(chrom.lengths) <- chrom.lengths.df$UCSC_seqlevel
  	} else {
  		names(chrom.lengths) <- chrom.lengths.df$NCBI_seqlevel
  	}
  	seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
  	
  	## Get counts at each position
  	ptm <- startTimedMessage("Getting counts ...")
  	data$methylated <- as.logical(as.integer(data$methylated)-1)
  	data$mask.dup <- c(FALSE, as.logical(data@seqnames[-1] == data@seqnames[-length(data)])) & c(FALSE, as.logical(data@ranges@start[-1] == data@ranges@start[-length(data)]))
  	data$unmeth.sum <- rev(cumsum(rev(!data$methylated)))
  	data$meth.sum <- rev(cumsum(rev(data$methylated)))
  	data.dedup <- data[!data$mask.dup]
  	data.dedup$unmeth.counts <- -diff(c(data.dedup$unmeth.sum, 0))
  	data.dedup$meth.counts <- -diff(c(data.dedup$meth.sum, 0))
  	stopTimedMessage(ptm)
  	
  	data.dedup$methylated <- NULL
  	data.dedup$mask.dup <- NULL
  	data.dedup$unmeth.sum <- NULL
  	data.dedup$meth.sum <- NULL
  	
  	## Add ratio and distance
  	data.dedup$ratio <- data.dedup$meth.counts / (data.dedup$unmeth.counts + data.dedup$meth.counts)
    data.dedup$distance <- c(-1, start(data.dedup)[-1] - end(data.dedup)[-length(data.dedup)] - 1)
    data.dedup$distance[data.dedup$distance < 0] <- Inf 
  	return(data.dedup)
}
