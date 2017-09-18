#' Methimpute data import
#' 
#' This page provides an overview of all \pkg{\link{methimpute}} data import functions.
#' @param file The file to import.
#' @param chrom.lengths A data.frame with chromosome names in the first, and chromosome lengths in the second column. Only chromosomes named in here will be returned. Alternatively a tab-separated file with such a data.frame (with headers).
#' @return A \code{\link{methimputeData}} object.
#' @name import
#' @examples 
#'## Get an example file in BSSeeker format
#'file <- system.file("extdata","arabidopsis_bsseeker.txt.gz", package="methimpute")
#'data(arabidopsis_chromosomes)
#'bsseeker.data <- importBSSeeker(file, chrom.lengths=arabidopsis_chromosomes)
#'
#'## Get an example file in Bismark format
#'file <- system.file("extdata","arabidopsis_bismark.txt", package="methimpute")
#'data(arabidopsis_chromosomes)
#'arabidopsis_chromosomes$chromosome <- sub('chr', '', arabidopsis_chromosomes$chromosome)
#'bismark.data <- importBismark(file, chrom.lengths=arabidopsis_chromosomes)
#'
NULL


#' @describeIn import Import a Methylpy methylation extractor file.
#' @importFrom utils read.table
#' @export
#' 
importMethylpy <- function(file, chrom.lengths=NULL) {
    
    ## Contexts
    contexts <- c(CG='CG.', CHG='C[ATC]G', CHH='C[ATC][ATC]')
    
    ## Import data
    classes <- c('character', 'numeric', 'character', 'character', 'numeric', 'numeric', 'numeric')
    ptm <- startTimedMessage("Reading file ", file, " ...")
    data.raw <- utils::read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
    data <- GRanges(seqnames=data.raw$V1, ranges=IRanges(start=data.raw$V2, end=data.raw$V2), strand=data.raw$V3, context=factor(NA, levels=names(contexts)), context.full=data.raw$V4)
    counts <- array(NA, dim=c(length(data), 2), dimnames=list(NULL, c("methylated", "total")))
    counts[,"methylated"] <- data.raw$V5
    counts[,"total"] <- data.raw$V6
    data$counts <- counts
    data$binomial.test <- data.raw$V7
    rm(data.raw)
    stopTimedMessage(ptm)
    
    ## Rework contexts
    ptm <- startTimedMessage("Reworking contexts ...")
    for (i1 in 1:length(contexts)) {
        ind <- grep(pattern = contexts[i1], x = data$context.full)
        data$context[ind] <- names(contexts)[i1]
    }
    data$context.full <- NULL
    stopTimedMessage(ptm)
    
    
    ## Assign seqlengths
    if (!is.null(chrom.lengths)) {
        if (is.character(chrom.lengths)) {
            df <- utils::read.table(chrom.lengths, header=TRUE)
        } else if (is.data.frame(chrom.lengths)) {
            df <- chrom.lengths
        }
        chrom.lengths <- df[,2]
        names(chrom.lengths) <- df[,1]
        # Filter by chromosomes supplied in chrom.lengths
        data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
        seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
    }
    
    return(data)
}


#' @describeIn import Import a BSSeeker methylation extractor file.
#' @importFrom utils read.table
#' @export
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
        } else if (is.data.frame(chrom.lengths)) {
            df <- chrom.lengths
        }
        chrom.lengths <- df[,2]
        names(chrom.lengths) <- df[,1]
        # Filter by chromosomes supplied in chrom.lengths
        data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
        seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
    }
    
    ## Make factors
    data$context <- factor(data$context)
    
    return(data)
}


#' @describeIn import Import a Bismark methylation extractor file.
#' @importFrom utils read.table
#' @export
#'
importBismark <- function(file, chrom.lengths=NULL) {
    
    # classes <- c(seqnames='character', position='numeric', strand='character', counts.methylated='numeric', counts.total='numeric', context='character', context.trinucleotide='character')
    classes <- c('character', 'numeric', 'character', 'numeric', 'numeric', 'character', 'character')
    ptm <- startTimedMessage("Reading file ", file, " ...")
    data.raw <- utils::read.table(file, skip=0, sep='\t', comment.char='', colClasses=classes)
    data <- GRanges(seqnames=data.raw$V1, ranges=IRanges(start=data.raw$V2, end=data.raw$V2), strand=data.raw$V3, context=data.raw$V6)
    counts <- array(NA, dim=c(length(data), 2), dimnames=list(NULL, c("methylated", "total")))
    counts[,"methylated"] <- data.raw$V4
    counts[,"total"] <- data.raw$V4 + data.raw$V5
    data$counts <- counts
    rm(data.raw)
    stopTimedMessage(ptm)
    
    
    ## Assign seqlengths
    if (!is.null(chrom.lengths)) {
        if (is.character(chrom.lengths)) {
            df <- utils::read.table(chrom.lengths, header=TRUE)
        } else if (is.data.frame(chrom.lengths)) {
            df <- chrom.lengths
        }
        chrom.lengths <- df[,2]
        names(chrom.lengths) <- df[,1]
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
#' @param chrom.lengths A data.frame with chromosome names in the first, and chromosome lengths in the second column. Only chromosomes named in here will be returned. Alternatively a tab-separated file with such a data.frame (with headers).
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
    
    seqlevels(data) <- paste0('chr', seqlevels(data))
    ## Assign seqlengths
    if (!is.null(chrom.lengths)) {
        if (is.character(chrom.lengths)) {
            df <- utils::read.table(chrom.lengths, header=TRUE)
        } else if (is.data.frame(chrom.lengths)) {
            df <- chrom.lengths
        }
        chrom.lengths <- df[,2]
        names(chrom.lengths) <- df[,1]
        # Filter by chromosomes supplied in chrom.lengths
        data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
        seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
    }
    
    ## Make factors
    data$context <- factor(data$context)
    
    return(data)
}
  
#' #' Import a Bismark methylation extractor file
#' #' 
#' #' Import a Bismark methylation extractor file into a \code{\link[GenomicRanges]{GRanges}} object.
#' #' 
#' #' @param files The files to import.
#' #' @param chrom.lengths A data.frame with chromosome names in the first, and chromosome lengths in the second column. Only chromosomes named in here will be returned. Alternatively a tab-separated file with such a data.frame (with headers).
#' #' @param temp.store A folder where to save temporary files on disk. Set to \code{NULL} if no temporary files should be created.
#' #' @return A \code{\link{methimputeData}} object.
#' #' 
#' #' @importFrom utils read.table
#' #' @export
#' importBismark <- function(files, chrom.lengths=NULL, temp.store=tempfile("importBismark")) {
#'   
#'     classes <- c('NULL', 'character', 'character', 'numeric', 'character')
#'     datas <- GRangesList()
#'     
#'     if (!is.null(temp.store)) {
#'         ## Create folder for temporary files
#'         if (!file.exists(temp.store)) {
#'             dir.create(temp.store)
#'         }
#'     }
#'     
#'     ### Do each file at a time to avoid memory overflow ###
#'     data.dedups <- GRangesList()
#'     for (i1 in 1:length(files)) {
#'         file <- files[i1]
#'         savename <- file.path(temp.store, paste0(basename(file), '.RData'))
#'         if (!file.exists(savename)) {
#'             ptm <- startTimedMessage("Reading file ", file, " ...")
#'             data.raw <- utils::read.table(file, skip=1, sep='\t', comment.char='', colClasses=classes)
#'             data <- GRanges(seqnames=data.raw$V3, ranges=IRanges(start=data.raw$V4, end=data.raw$V4), strand="*", methylated=factor(data.raw$V2, levels=c("-", "+")), context=data.raw$V5)
#'             rm(data.raw)
#'             stopTimedMessage(ptm)
#'             if (!is.null(temp.store)) {
#'                 ptm <- startTimedMessage("Saving to file ", savename, " ...")
#'                 save(data, file=savename)
#'                 stopTimedMessage(ptm)
#'             }
#'         } else {
#'             ptm <- startTimedMessage("Loading file ", savename, " ...")
#'             temp.env <- new.env()
#'             data <- get(load(savename, envir=temp.env), envir=temp.env)
#'             stopTimedMessage(ptm)
#'         }
#'         
#'         ## Sorting
#'         ptm <- startTimedMessage("Sorting ...")
#'         data <- sort(data)
#'         stopTimedMessage(ptm)
#'         
#'         ## Get counts at each position
#'         ptm <- startTimedMessage("Getting counts ...")
#'         data$methylated <- as.logical(as.integer(data$methylated)-1)
#'         data$mask.dup <- c(FALSE, as.logical(data@seqnames[-1] == data@seqnames[-length(data)])) & c(FALSE, as.logical(data@ranges@start[-1] == data@ranges@start[-length(data)]))
#'         data$unmeth.sum <- rev(cumsum(rev(!data$methylated)))
#'         data$meth.sum <- rev(cumsum(rev(data$methylated)))
#'         data.dedup <- data[!data$mask.dup]
#'         rm(data)
#'         data.dedup$counts.unmethylated <- -diff(c(data.dedup$unmeth.sum, 0))
#'         data.dedup$counts.methylated <- -diff(c(data.dedup$meth.sum, 0))
#'         data.dedup$methylated <- NULL
#'         data.dedup$mask.dup <- NULL
#'         data.dedup$unmeth.sum <- NULL
#'         data.dedup$meth.sum <- NULL
#'         data.dedups[[length(data.dedups)+1]] <- data.dedup
#'         stopTimedMessage(ptm)
#'         
#'         ## Merging counts from data objects
#'         if (length(data.dedups) == 2) {
#'             ptm <- startTimedMessage("Merging counts ...")
#'             data <- unlist(data.dedups, use.names=FALSE)
#'             data <- sort(data)
#'             rm(data.dedups)
#'             data$mask.dup <- c(FALSE, as.logical(data@seqnames[-1] == data@seqnames[-length(data)])) & c(FALSE, as.logical(data@ranges@start[-1] == data@ranges@start[-length(data)]))
#'             data$unmeth.sum <- rev(cumsum(rev(data$counts.unmethylated)))
#'             data$meth.sum <- rev(cumsum(rev(data$counts.methylated)))
#'             data.dedup <- data[!data$mask.dup]
#'             rm(data)
#'             data.dedup$counts.unmethylated <- -diff(c(data.dedup$unmeth.sum, 0))
#'             data.dedup$counts.methylated <- -diff(c(data.dedup$meth.sum, 0))
#'             data.dedup$mask.dup <- NULL
#'             data.dedup$unmeth.sum <- NULL
#'             data.dedup$meth.sum <- NULL
#'             data.dedups <- GRangesList(data.dedup)
#'             stopTimedMessage(ptm)
#'         }
#'         
#'     }
#'     data <- data.dedups[[1]]
#'     rm(data.dedups, data.dedup)
#'   
#'     ## Assign seqlengths
#'     if (!is.null(chrom.lengths)) {
#'         if (is.character(chrom.lengths)) {
#'             df <- utils::read.table(chrom.lengths, header=TRUE)
#'         } else if (is.data.frame(chrom.lengths)) {
#'             df <- chrom.lengths
#'         }
#'         chrom.lengths <- df[,2]
#'         names(chrom.lengths) <- df[,1]
#'         # Filter by chromosomes supplied in chrom.lengths
#'         data <- keepSeqlevels(data, seqlevels(data)[seqlevels(data) %in% names(chrom.lengths)])
#'         seqlengths(data) <- chrom.lengths[names(seqlengths(data))]
#'     }
#'     
#'     ## Make factors
#'     data$context <- factor(data$context)
#'     
#'     return(data)
#' }
#' 
