exportGRangesAsBed <- function(gr, grlist=NULL, filename, namecol=NULL, scorecol=NULL, colorcol=NULL, trackname=filename, description=trackname, colors=NULL, priority=0, header=TRUE) {
  
    ## Concatenate into list ##
    if (is.null(grlist)) {
        grlist <- list(gr)
        names(grlist)[1] <- trackname
    }
    if (is.null(names(grlist))) {
        stop("'grlist' must be named")
    }
  
    ## Generate colors for export ##
    itemRGB.string <- 'Off'
    if (!is.null(colorcol)) {
        itemRGB.string <- 'On'
        unique.names <- list()
        for (i1 in 1:length(grlist)) {
            unique.names[[i1]] <- levels(mcols(grlist[[i1]])[,colorcol])
            if (is.null(unique.names[[i1]])) {
                unique.names[[i1]] <- unique(mcols(grlist[[i1]])[,colorcol])
            }
        }
        unique.names <- unlist(unique.names)
        unique.names <- unique.names[!duplicated(unique.names)]
        color.mapping <- getDistinctColors(unique.names)
        if (length(color.mapping) == 1) {
            track.color.mapping <- paste(grDevices::col2rgb(color.mapping), collapse=',')
        } else if (length(color.mapping) > 1) {
            track.color.mapping <- apply(grDevices::col2rgb(color.mapping), 2, paste, collapse=',')
        }
    } else if (!is.null(colors)) {
        itemRGB.string <- 'On'
        if (length(colors) == 1) {
            track.colors <- paste(grDevices::col2rgb(colors), collapse=',')
        } else {
            track.colors <- apply(grDevices::col2rgb(colors), 2, paste, collapse=',')
        }
        track.colors <- rep(track.colors, length=length(grlist))[1:length(grlist)]
    }
  
    ## Open file to write to ##
    filename.bed.gz <- paste0(filename,".bed.gz")
    filehandle <- gzfile(filename.bed.gz, 'w')
    ptm <- startTimedMessage("Writing to file ", filename.bed.gz, " ...")
    cat("", file=filehandle)
  
    ## Write tracks to file ##
    for (i1 in 1:length(grlist)) {
        trackname.string <- paste0(names(grlist)[i1], ': ', trackname)
        description.string <- paste0(names(grlist)[i1], ': ', description)
        gri <- grlist[[i1]]
        gri <- insertchr(gri)
        # Convert from 1-based closed to 0-based half open
        df <- data.frame(seqnames=gri$chromosome, start=start(gri)-1, end=end(gri), name=".", score=1000, strand=sub('\\*', '.', strand(gri)))
        if (!is.null(namecol)) {
            df$name <- mcols(gri)[,namecol]
        }
        if (!is.null(scorecol)) {
            df$score <- mcols(gri)[,scorecol]
        }
        if (!is.null(colorcol)) {
            df$thickStart <- df$start
            df$thickEnd <- df$end
            df$itemRgb <- track.color.mapping[as.character(mcols(gri)[,colorcol])]
        } else if (!is.null(colors)) {
            df$thickStart <- df$start
            df$thickEnd <- df$end
            df$itemRgb <- track.colors[i1]
        }
        if (header) {
            cat(paste0('track name="', trackname.string,'" description="', description.string,'" visibility=1 itemRgb=', itemRGB.string, ' priority=', priority, '\n'), file=filehandle, append=TRUE)
        }
        if (nrow(df) == 0) {
            warning("Nothing to export in variable 'grlist[[", i1, "]]'.")
        } else {
            utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filehandle, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        }
    }
    
    ## Close file ##
    close(filehandle)
    stopTimedMessage(ptm)
  
}



exportGRangesAsWiggle <- function(gr, grlist=NULL, col2export, filename, trackname=filename, description=trackname, colors='gray68', priority=0) {
  
    ## Concatenate into list ##
    if (is.null(grlist)) {
        grlist <- list(gr)
        names(grlist)[1] <- trackname
    }
    if (is.null(names(grlist))) {
        stop("'grlist' must be named")
    }
  
    ## Generate colors for export ##
    if (length(colors) == 1) {
        track.colors <- paste(grDevices::col2rgb(colors), collapse=',')
    } else {
        track.colors <- apply(grDevices::col2rgb(colors), 2, paste, collapse=',')
    }
    track.colors <- rep(track.colors, length=length(grlist))[1:length(grlist)]
    
    ## Open file to write to ##
    filename.wig.gz <- paste0(filename,".wig.gz")
    filehandle <- gzfile(filename.wig.gz, 'w')
    ptm <- startTimedMessage("Writing to file ", filename.wig.gz, " ...")
    cat("", file=filehandle)
    
    ## Write tracks to file ##
    for (i1 in 1:length(grlist)) {
        trackname.string <- paste0(names(grlist)[i1], ': ', trackname)
        description.string <- paste0(names(grlist)[i1], ': ', description)
        gri <- grlist[[i1]]
        gri <- insertchr(gri)
        cat(paste0('track type=wiggle_0 name="', trackname.string,'" description="', description.string,'" visibility=full autoScale=on color=', track.colors[i1],' maxHeightPixels=100:50:20 graphType=bar priority=', priority,'\n'), file=filehandle, append=TRUE)
        for (chrom in unique(gri$chromosome)) {
            gri.chrom <- gri[gri$chromosome == chrom]
            binsize <- width(gri.chrom)[1]
            cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filehandle, append=TRUE)
            utils::write.table(mcols(gri.chrom)[,col2export], file=filehandle, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
    }
    
    ## Close file ##
    close(filehandle)
    stopTimedMessage(ptm)
  
}

  
#=====================================================
# Helper functions   
#=====================================================
insertchr <- function(gr) {
    # Change chromosome names from '1' to 'chr1' if necessary
    mask <- which(!grepl('chr', seqnames(gr)))
    mcols(gr)$chromosome <- as.character(seqnames(gr))
    mcols(gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(gr)$chromosome[mask])
    mcols(gr)$chromosome <- as.factor(mcols(gr)$chromosome)
    return(gr)
}

