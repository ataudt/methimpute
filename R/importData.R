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

#' Import a BSSeeker methylation extractor file
#' 
#' Import a BSSeeker methylation extractor file into a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param file The file to import.
#' @param chrom.lengths A named vector containing the chromosome lengths. Only chromosomes named in here will be returned.
#' @return A \code{\link{methimputeData}} object.
#' 
#' @importFrom utils read.table
#' @export
#' 
#' @examples 
#'## Get an example file in BSSeeker format
#'file <- system.file("extdata","arabidopsis_bsseeker.txt.gz", package="methimpute")
#'data(arabidopsis_chromosomes)
#'bsseeker.data <- importBSSeeker(file, chrom.lengths=arabidopsis_chromosomes)
#'
importBSSeeker <- function(file, chrom.lengths=NULL) {
    
    # classes <- c(seqnames='character', nucleotide='character', position='numeric', context='character', context.dinucleotide='character', methylation.level='numeric', counts.methylated='numeric', counts.total='numeric')
    classes <- c('character', 'character', 'numeric', 'character', 'character', 'numeric', 'numeric', 'numeric')
    ptm <- startTimedMessage("Reading file ", file, " ...")
	  data.raw <- utils::read.table(file, skip=0, sep='\t', comment.char='', colClasses=classes)
	  data <- GRanges(seqnames=data.raw$V1, ranges=IRanges(start=data.raw$V3, end=data.raw$V3), strand=c('C'='+', 'G'='-')[data.raw$V2], context=data.raw$V4)
  	counts <- array(NA, dim=c(length(data), 2), dimnames=list(NULL, c("methylated", "total")))
  	counts[,"methylated"] <- data.raw$V7
  	counts[,"total"] <- data.raw$V8
  	data$counts <- counts
  	rm(data.raw)
  	stopTimedMessage(ptm)
	  
  	
  	## Assign seqlengths
  	if (!is.null(chrom.lengths)) {
  	    if (is.character(chrom.lengths)) {
  	        df <- utils::read.table(chrom.lengths, header=TRUE)
  	        chrom.lengths <- df[,2]
  	        names(chrom.lengths) <- df[,1]
  	    }
      	# Filter by chromosomes supplied in chrom.lengths
      	data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
      	seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
  	}
  	
    ## Make factors
    data$context <- factor(data$context)
    
  	return(data)
}


#' Import a Rene methylation extractor file
#' 
#' Import a Rene methylation extractor file into a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param file The file to import.
#' @param chrom.lengths A named vector containing the chromosome lengths. Only chromosomes named in here will be returned.
#' @return A \code{\link{methimputeData}} object.
#' 
#' @importFrom utils read.table
importRene <- function(file, chrom.lengths=NULL) {
  
  	classes <- c('character', 'numeric', 'character', 'NULL', 'character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric')
    ptm <- startTimedMessage("Reading file ", file, " ...")
	  data.raw <- utils::read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
  	data <- GRanges(seqnames=data.raw$V1, ranges=IRanges(start=data.raw$V2, end=data.raw$V2), strand=c('F'='+', 'R'='-')[data.raw$V3], methylated=data.raw$V10, context=data.raw$V5)
  	counts <- array(NA, dim=c(length(data), 2), dimnames=list(NULL, c("methylated", "total")))
  	counts[,"methylated"] <- data.raw$V6
  	counts[,"total"] <- data.raw$V7
  	data$counts <- counts
  	rm(data.raw)
  	stopTimedMessage(ptm)
  	
  	## Assign seqlengths
  	if (!is.null(chrom.lengths)) {
  	    if (is.character(chrom.lengths)) {
  	        df <- utils::read.table(chrom.lengths, header=TRUE)
  	        chrom.lengths <- df[,2]
  	        names(chrom.lengths) <- df[,1]
  	    }
      	# Filter by chromosomes supplied in chrom.lengths
      	data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
      	seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
  	}
  	
    ## Make factors
    data$context <- factor(data$context)
    
  	return(data)
}
  
#' Import a Bismarck methylation extractor file
#' 
#' Import a Bismarck methylation extractor file into a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param files The files to import.
#' @param chrom.lengths A named vector containing the chromosome lengths. Only chromosomes named in here will be returned.
#' @param temp.store A folder where to save temporary files on disk. Set to \code{NULL} if no temporary files should be created.
#' @return A \code{\link{methimputeData}} object.
#' 
#' @importFrom utils read.table
#' @export
importBismarck <- function(files, chrom.lengths=NULL, temp.store=tempfile("importBismarck")) {
  
  	classes <- c('NULL', 'character', 'character', 'numeric', 'character')
  	datas <- GRangesList()
  	
	  if (!is.null(temp.store)) {
        ## Create folder for temporary files
      	if (!file.exists(temp.store)) {
      	    dir.create(temp.store)
      	}
	  }
  	
  	### Do each file at a time to avoid memory overflow ###
  	data.dedups <- GRangesList()
  	for (i1 in 1:length(files)) {
  	    file <- files[i1]
	      savename <- file.path(temp.store, paste0(basename(file), '.RData'))
	      if (!file.exists(savename)) {
      	    ptm <- startTimedMessage("Reading file ", file, " ...")
        	  data.raw <- utils::read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
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
      	data.dedup$counts.unmethylated <- -diff(c(data.dedup$unmeth.sum, 0))
      	data.dedup$counts.methylated <- -diff(c(data.dedup$meth.sum, 0))
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
          	data$unmeth.sum <- rev(cumsum(rev(data$counts.unmethylated)))
          	data$meth.sum <- rev(cumsum(rev(data$counts.methylated)))
          	data.dedup <- data[!data$mask.dup]
          	rm(data)
          	data.dedup$counts.unmethylated <- -diff(c(data.dedup$unmeth.sum, 0))
          	data.dedup$counts.methylated <- -diff(c(data.dedup$meth.sum, 0))
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
  	if (!is.null(chrom.lengths)) {
  	    if (is.character(chrom.lengths)) {
  	        df <- utils::read.table(chrom.lengths, header=TRUE)
  	        chrom.lengths <- df[,2]
  	        names(chrom.lengths) <- df[,1]
  	    }
      	# Filter by chromosomes supplied in chrom.lengths
      	data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
      	seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
  	}
  	
    ## Make factors
    data$context <- factor(data$context)
    
  	return(data)
}

