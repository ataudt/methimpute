#' Extract cytosine coordinates
#' 
#' Extract cytosine coordinates and context information from a FASTA file. Cytosines in ambiguous contexts are not reported.
#' 
#' @param file A character with the file name.
#' @param contexts The contexts that should be extracted.
#' 
#' @importFrom Biostrings readDNAStringSet vmatchPattern reverseComplement
#' @export
extractCytosinesFromFASTA <- function(file, contexts = c('CG','CHG','CHH')) {
  
    ### Read file
    fasta <- Biostrings::readDNAStringSet(file)
    
    ### Chromosome lengths
    chromlengths <- width(fasta)
    names(chromlengths) <- names(fasta)
    
    ### Positions with ambiguous nucleotide codes
    ptm <- startTimedMessage("Scanning for ambiguous nucleotides ...")
    notuse <- GRangesList()
    ambicodes <- c('R','Y','S','W','K','M','B','D','H','V','N')
    for (ambicode in ambicodes) {
        notuse.string <- Biostrings::vmatchPattern(ambicode, subject = fasta, fixed = TRUE)
        notuse[[ambicode]] <- as(notuse.string, 'GRanges')
    }
    notuse <- unlist(notuse, use.names = FALSE)
    stopTimedMessage(ptm)
    
    ### Extract cytosine positions with context
    ptm <- startTimedMessage("Extracting cytosines from forward strand ...")
    cytosines.forward <- GRangesList()
    for (context in contexts) {
        positions <- Biostrings::vmatchPattern(context, subject = fasta, fixed = FALSE)
        positions <- as(positions, 'GRanges')
        strand(positions) <- '+'
        positions$context <- factor(context, levels=contexts)
        cytosines.forward[[context]] <- positions
    }
    cytosines.forward <- unlist(cytosines.forward, use.names = FALSE)
    seqlengths(cytosines.forward) <- chromlengths
    # Remove ambiguous positions
    ind <- findOverlaps(cytosines.forward, notuse)
    cytosines.forward <- cytosines.forward[-ind@from]
    # Set width to 1
    end(cytosines.forward) <- start(cytosines.forward)
    stopTimedMessage(ptm)
    
    ### Extract reverse context
    ptm <- startTimedMessage("Extracting cytosines from reverse strand ...")
    fasta <- Biostrings::reverseComplement(fasta)
    cytosines.reverse <- GRangesList()
    for (context in contexts) {
        positions <- Biostrings::vmatchPattern(context, subject = fasta, fixed = FALSE)
        positions <- as(positions, 'GRanges')
        strand(positions) <- '-'
        positions$context <- factor(context, levels=contexts)
        cytosines.reverse[[context]] <- positions
    }
    cytosines.reverse <- unlist(cytosines.reverse, use.names = FALSE)
    seqlengths(cytosines.reverse) <- chromlengths
    # Correct coordinates from reverseComplement
    starts <- chromlengths[as.character(seqnames(cytosines.reverse))] - start(cytosines.reverse) + 1
    ends <- chromlengths[as.character(seqnames(cytosines.reverse))] - end(cytosines.reverse) + 1
    cytosines.reverse <- GRanges(seqnames = seqnames(cytosines.reverse), ranges = IRanges(start = ends, end = starts), strand = '-', context = cytosines.reverse$context)
    # Remove ambiguous positions
    ind <- findOverlaps(cytosines.reverse, notuse)
    cytosines.reverse <- cytosines.reverse[-ind@from]
    # Set width to 1
    start(cytosines.reverse) <- end(cytosines.reverse)
    stopTimedMessage(ptm)
    
    ### Merge and sort
    ptm <- startTimedMessage("Merging ...")
    cytosines <- c(cytosines.forward, cytosines.reverse)
    cytosines <- sort(cytosines, ignore.strand=TRUE)
    stopTimedMessage(ptm)
    
    return(cytosines)
}