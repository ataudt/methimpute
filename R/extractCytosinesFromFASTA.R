#' Extract cytosine coordinates
#' 
#' Extract cytosine coordinates and context information from a FASTA file. Cytosines in ambiguous reference contexts are not reported.
#' 
#' @param file A character with the file name.
#' @param contexts The contexts that should be extracted.
#' @param anchor.C A named vector with positions of the anchoring C in the \code{contexts}. This is necessary to distinguish contexts such as C*C*CG (anchor.C = 2) and *C*CCG (anchor.C = 1). Names must match the contexts. If unspecified, the first C within each context will be taken as anchor.
#' 
#' @importFrom Biostrings readDNAStringSet vmatchPattern reverseComplement
#' @export
#' @examples
#' ## Read a non-compressed FASTA files:
#' filepath <- system.file("extdata", "arabidopsis_sequence.fa.gz", package="methimpute")
#' cytosines <- extractCytosinesFromFASTA(filepath)
#' 
#' ## Split CG context into subcontexts
#' cytosines <- extractCytosinesFromFASTA(filepath, contexts = 'CG')
#' cytosines <- extractCytosinesFromFASTA(filepath,
#'                contexts = c('DCG', 'CCG'),
#'                anchor.C = c(DCG=2, CCG=2))
extractCytosinesFromFASTA <- function(file, contexts = c('CG','CHG','CHH'), anchor.C = NULL) {
  
    ### C anchors ###
    if (is.null(anchor.C)) {
        anchor.C <- regexpr('C', contexts)[1:length(contexts)]
        names(anchor.C) <- contexts
    }
    if (!all(names(anchor.C)==contexts)) {
        stop("names(anchor.C) must be equal to contexts")
    }
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
    if (length(ind) > 0) {
        cytosines.forward <- cytosines.forward[-ind@from]
    }
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
    if (length(ind) > 0) {
        cytosines.reverse <- cytosines.reverse[-ind@from]
    }
    # Set width to 1
    start(cytosines.reverse) <- end(cytosines.reverse)
    stopTimedMessage(ptm)
    
    ### Merge and sort
    ptm <- startTimedMessage("Merging ...")
    cytosines <- c(cytosines.forward, cytosines.reverse)
    cytosines <- sort(cytosines, ignore.strand=TRUE)
    stopTimedMessage(ptm)
    
    # Shift positions by position of anchor C
    cind <- anchor.C[cytosines$context]
    strandint <- c('-'=-1,'+'=1,'*'=1)[as.character(strand(cytosines))]
    starts <- start(cytosines) + strandint * cind - strandint * 1
    ends <- end(cytosines) + strandint * cind - strandint * 1
    cytosines <- GRanges(seqnames=seqnames(cytosines), ranges=IRanges(start=starts, end=ends), strand=strand(cytosines), context=cytosines$context)
    
    return(cytosines)
}