#' Enrichment analysis
#' 
#' Plotting functions for enrichment analysis of \code{\link{multimodel}} or \code{\link{combinedMultimodel}} objects with any annotation of interest, specified as a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @name enrichment_analysis
#' @param combinations A vector with combinations for which the enrichment will be calculated. If \code{NULL} all combinations will be considered.
#' @param marks A vector with marks for which the enrichment is plotted. If \code{NULL} all marks will be considered.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object containing the plot or a list() with \code{\link[ggplot2:ggplot]{ggplot}} objects if several plots are returned. For \code{plotFoldEnrichHeatmap} a named array with fold enrichments if \code{plot=FALSE}.
#' @author Aaron Taudt
#' @seealso \code{\link{plotting}}
#' @examples 
#'### Get an example multimodel ###
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'
#'### Obtain gene coordinates for rat from biomaRt ###
#'library(biomaRt)
#'ensembl <- useMart('ENSEMBL_MART_ENSEMBL', host='may2012.archive.ensembl.org',
#'                   dataset='rnorvegicus_gene_ensembl')
#'genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position',
#'                            'end_position', 'strand', 'external_gene_id',
#'                            'gene_biotype'),
#'               mart=ensembl)
#'# Transform to GRanges for easier handling
#'genes <- GRanges(seqnames=paste0('chr',genes$chromosome_name),
#'                 ranges=IRanges(start=genes$start, end=genes$end),
#'                 strand=genes$strand,
#'                 name=genes$external_gene_id, biotype=genes$gene_biotype)
#'print(genes)
#'
#'### Make the enrichment plots ###
#'# We expect promoter [H3K4me3] and bivalent-promoter signatures [H3K4me3+H3K27me3]
#'# to be enriched at transcription start sites.
#'    plotEnrichment(model = model, annotation = genes, bp.around.annotation = 15000) +
#'    ggtitle('Fold enrichment around genes') +
#'    xlab('distance from gene body')
#'  
#'# Plot enrichment only at TSS. We make use of the fact that TSS is the start of a gene.
#'    plotEnrichment(model, genes, region = 'start') +
#'    ggtitle('Fold enrichment around TSS') +
#'    xlab('distance from TSS in [bp]')
#'# Note: If you want to facet the plot because you have many combinatorial states you
#'# can do that with
#'    plotEnrichment(model, genes, region = 'start') +
#'    facet_wrap(~ combination)
#'  
#'# Another form of visualization that shows every TSS in a heatmap
#'# If transparency is not supported try to plot to pdf() instead.
#'    tss <- resize(genes, width = 3, fix = 'start')
#'    plotEnrichCountHeatmap(model, tss) +
#'    theme(strip.text.x = element_text(size=6))
#'  
#'# Fold enrichment with different biotypes, showing that protein coding genes are
#'# enriched with (bivalent) promoter combinations [H3K4me3] and [H3K4me3+H3K27me3],
#'# while rRNA is enriched with the empty [] and repressive combinations [H3K27me3].
#'    tss <- resize(genes, width = 3, fix = 'start')
#'    biotypes <- split(tss, tss$biotype)
#'    plotFoldEnrichHeatmap(model, annotations=biotypes) + coord_flip()
#'
NULL


#' @describeIn enrichment_analysis Compute the fold enrichment of combinatorial states for multiple annotations.
#' @param model A \code{\link{combinedMultimodel}} or \code{\link{multimodel}} object or a file that contains such an object.
#' @param annotations A \code{list()} with \code{\link{GRanges}} objects containing coordinates of multiple annotations The names of the list entries will be used to name the return values.
#' @param plot A logical indicating whether the plot or an array with the fold enrichment values is returned.
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges subsetByOverlaps
#' @importFrom reshape2 melt
#' @export
heatmapFoldEnrichment <- function(model, annotations, plot=TRUE) {
    
    ## Variables
    bins <- model$data
    genome <- sum(as.numeric(width(bins)))
    annotationsAtBins <- lapply(annotations, function(x) { IRanges::subsetByOverlaps(x, bins) })
    feature.lengths <- lapply(annotationsAtBins, function(x) { sum(as.numeric(width(x))) })
    state.levels <- levels(bins$state)
    # Only state levels that actually appear
    state.levels <- state.levels[state.levels %in% unique(bins$state)]
    
    ggplts <- list()
    folds <- list()
    maxfolds <- list()
    
    ptm <- startTimedMessage("Calculating fold enrichment ...")
    fold <- array(NA, dim=c(length(annotationsAtBins), length(state.levels)), dimnames=list(annotation=names(annotationsAtBins), state=state.levels))
    for (istate in 1:length(state.levels)) {
        mask <- bins$state == state.levels[istate]
        bins.mask <- bins[mask]
        state.length <- sum(as.numeric(width(bins.mask)))
        for (ifeat in 1:length(annotationsAtBins)) {
            feature <- annotationsAtBins[[ifeat]]
            ind <- findOverlaps(bins.mask, feature)

            binsinfeature <- bins.mask[unique(S4Vectors::queryHits(ind))]
            sum.binsinfeature <- sum(as.numeric(width(binsinfeature)))

            featuresinbins <- feature[unique(S4Vectors::subjectHits(ind))]
            sum.featuresinbins <- sum(as.numeric(width(featuresinbins)))

            fold[ifeat,istate] <- sum.binsinfeature / state.length / feature.lengths[[ifeat]] * genome
        }
    }
    stopTimedMessage(ptm)
    
    ## Transform fold values to range centered around 0
    logfold <- log(fold)
    
    if (plot) {
        logfold[is.infinite(logfold)] <- NaN
        colorbar <- gplots::colorpanel(n = 100, low="blue", mid="white",high="red")
        charf <- 0.6
        margins <- c( max(5, charf*max(nchar(colnames(logfold)))) , max(5, charf*max(nchar(rownames(logfold)))))
        gplots::heatmap.2(logfold, col=colorbar, trace='none', density.info='none', key.title = 'log(observed/expected)', key.xlab='log(observed/expected)', margins=margins, na.color='gray')
    } else {
        return(logfold)
    }
    
    # ## Clustering
    # if (cluster) {
    #   hc.anno <- hclust(dist(fold))
    #   hc.state <- hclust(dist(t(fold)))
    #   fold <- fold[hc.anno$order, hc.state$order]
    # }
    # 
    # if (plot) {
    #   ## Dendrograms
    #   if (cluster) {
    #     # Row dendrogram
    #     row.dendro <- ggplot(ggdendro::segment(ggdendro::dendro_data(as.dendrogram(hc.anno)))) +  geom_segment(aes_string(x='x', y='y', xend='xend', yend='yend')) +  theme_none + theme(axis.title.x=element_blank()) + coord_flip()
    #     # Column dendrogram
    #     col.dendro <- ggplot(ggdendro::segment(ggdendro::dendro_data(as.dendrogram(hc.state)))) +  geom_segment(aes_string(x='x', y='y', xend='xend', yend='yend')) + theme_none
    #   }
    #   df <- reshape2::melt(fold, value.name='foldEnrichment')
    #   df$state <- factor(df$state, levels=colnames(fold))
    #   maxfold <- max(df$foldEnrichment, na.rm=TRUE)
    #   ggplt <- ggplot(df) + geom_tile(aes_string(x='state', y='annotation', fill='foldEnrichment'))
    #   ggplt <- ggplt + theme_bw()
    #   ggplt <- ggplt + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
    #   ggplt <- ggplt + scale_fill_gradientn(colors=grDevices::colorRampPalette(c("blue","white","red"))(20), values=c(seq(0,1,length.out=10), seq(1,maxfold,length.out=10)), rescaler=function(x,...) {x}, oob=identity, limits=c(0,maxfold))
    #   cowplt <- cowplot::plot_grid(NULL, col.dendro, row.dendro, ggplt, align = 'hv', ncol = 2)
    #   return(ggplt)
    # } else {
    #   return(fold)
    # }
}


#' @describeIn enrichment_analysis Plot read counts around annotation as heatmap.
#' @inheritParams enrichmentAtAnnotation
#' @param max.rows An integer specifying the number of randomly subsampled rows that are plotted from the \code{annotation} object. This is necessary to avoid crashing for heatmaps with too many rows.
#' @importFrom reshape2 melt
#' @importFrom IRanges subsetByOverlaps
#' @export
plotEnrichCountHeatmap <- function(model, annotation, bp.around.annotation=10000, max.rows=1000) {

    ## Variables
    bins <- model$bins
    if (class(model) == class.combined.multivariate.model) {
        conditions <- sub('combination.', '', grep('combination', names(mcols(bins)), value=TRUE))
        comb.levels <- levels(mcols(bins)[,paste0('combination.', conditions[1])])
        ## Create new column combination with all conditions combined
        combs <- list()
        for (condition in conditions) {
            combs[[condition]] <- paste0(condition, ":", mcols(bins)[,paste0('combination.', condition)])
        }
        combs$sep <- ', '
        bins$combination <- factor(do.call(paste, combs))
    } else if (class(model) == class.multivariate.model) {
        comb.levels <- levels(bins$combination)
    }
    binsize <- width(bins)[1]
    around <- round(bp.around.annotation/binsize)
    
    ## Get RPKM values
    bins$counts <- sweep(bins$counts, MARGIN = 2, STATS = colSums(bins$counts), FUN = '/')
    bins$counts <- bins$counts * 1e6 * 1000/mean(width(bins))

    # Subsampling for plotting of huge data.frames
    annotation <- IRanges::subsetByOverlaps(annotation, bins)
    if (length(annotation)>max.rows) {
        annotation <- sample(annotation, size=max.rows, replace=FALSE)
    }
  
    # Get bins that overlap the start of annotation
    ptm <- startTimedMessage("Overlaps with annotation ...")
    index.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], bins, select="first")
    index.minus <- findOverlaps(annotation[strand(annotation)=='-'], bins, select="last")
    index.plus <- index.plus[!is.na(index.plus)]
    index.minus <- index.minus[!is.na(index.minus)]
    index <- c(index.plus, index.minus)
    stopTimedMessage(ptm)
    
    # Get surrounding indices
    ptm <- startTimedMessage("Getting surrounding indices ...")
    ext.index.plus <- array(NA, dim=c(length(index.plus), 2*around+1), dimnames=list(anno=1:length(index.plus), position=binsize*seq(-around, around, 1)))
    for (i1 in 1:length(index.plus)) {
        ext.index.plus[i1,] <- seq(from=-around+index.plus[i1], to=index.plus[i1]+around)
    }
    if (length(index.minus)>0) {
        ext.index.minus <- array(NA, dim=c(length(index.minus), 2*around+1), dimnames=list(anno=1:length(index.minus), position=binsize*seq(-around, around, 1)))
        for (i1 in 1:length(index.minus)) {
            ext.index.minus[i1,] <- rev( seq(from=-around+index.minus[i1], to=index.minus[i1]+around) )
        }
        ext.index <- rbind(ext.index.plus, ext.index.minus)
    } else {
        ext.index <- ext.index.plus
    }
    ext.index[ext.index <= 0] <- NA
    ext.index[ext.index > length(bins)] <- NA
    rownames(ext.index) <- 1:nrow(ext.index)
    stopTimedMessage(ptm)
    
    ## Go through combinations and then tracks to get the read counts
    ptm <- startTimedMessage("Getting read counts")
    counts <- list()
    combinations <- names(sort(table(bins$combination[index]), decreasing = TRUE))
    for (combination in combinations) {
        counts[[combination]] <- list()
        index.combination <- which(bins$combination[index]==combination)
        ext.index.combination <- ext.index[index.combination,]
        if (is.null(dim(ext.index.combination))) {
            ext.index.combination <- array(ext.index.combination, dim=c(1,dim(ext.index)[[2]]), dimnames=list(anno=rownames(ext.index)[index.combination], position=dimnames(ext.index)[[2]]))
        }
        for (ntrack in colnames(bins$counts)) {
            counts[[combination]][[ntrack]] <- array(bins$counts[ext.index.combination,ntrack], dim=dim(ext.index.combination), dimnames=dimnames(ext.index.combination))
        }
    }
    stopTimedMessage(ptm)
    
    ## Theme
    custom_theme <- theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill='NA'),
        panel.background = element_rect(fill='white'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
    )
    
    ## Prepare data.frame
    ptm <- startTimedMessage("Making the plot ...")
    # Exclude rare combinations for plotting
    num.comb <- sapply(counts, function(x) { nrow(x[[1]]) })
    comb2keep <- names(num.comb)[num.comb/sum(num.comb) > 0.005]
    counts <- counts[comb2keep]
    df <- reshape2::melt(counts)
    names(df) <- c('id','position','RPKM','track','combination')
    df$id <- factor(df$id, levels=rev(unique(df$id)))
    df$combination <- factor(df$combination, levels=unique(df$combination))
    df$track <- factor(df$track, levels=colnames(bins$counts))
    
    ## Plot as heatmap
    ggplt <- ggplot(df) + geom_tile(aes_string(x='position', y='id', color='combination'))
    ggplt <- ggplt + scale_color_manual(values = getDistinctColors(length(unique(df$combination))))
    ggplt <- ggplt + geom_tile(aes_string(x='position', y='id', fill='RPKM'), alpha=0.6)
    ggplt <- ggplt + facet_wrap( ~ track, nrow=1) + custom_theme
    ggplt <- ggplt + xlab('distance from annotation in [bp]') + ylab('')
    ggplt <- ggplt + scale_fill_continuous(trans='log1p', low='white', high='black')
    # Insert horizontal lines
    y.lines <- sapply(split(df$id, df$combination), function(x) { max(as.integer(x)) })
    df.lines <- data.frame(y=sort(y.lines[-1]) + 0.5)
    ggplt <- ggplt + geom_hline(data=df.lines, mapping=aes_string(yintercept='y'), linetype=2)
    # Increase color size in legend
    ggplt <- ggplt + guides(color=guide_legend(override.aes = list(size=2)))
    stopTimedMessage(ptm)
    
    return(ggplt)
}


#' @describeIn enrichment_analysis Plot fold enrichment of combinatorial states around and inside of annotation.
#' @inheritParams enrichmentAtAnnotation
#' @importFrom reshape2 melt
#' @export
plotEnrichment <- function(model, annotation, bp.around.annotation=10000, region=c("start","inside","end"), num.intervals=20, what='combinations', combinations=NULL, marks=NULL, statistic='fold') {

    ## Check user input
    if ((!what %in% c('combinations','peaks','counts')) | length(what) > 1) {
        stop("argument 'what' must be one of c('combinations','peaks','counts')")
    }
    if (!is.null(marks) & what!='peaks') {
        stop("Please set argument 'what=\"peaks\"' if you want to plot marks instead of combinations.")
    }
  
    ## Variables
    model$bins$counts <- rpkm.matrix(model$bins$counts, binsize=mean(width(model$bins)))
    bins <- model$bins
    if (class(model) == class.combined.multivariate.model) {
    } else if (class(model) == class.multivariate.model) {
        # Rename 'combination' to 'combination.' for coherence with combinedMultimodel
        names(mcols(bins))[grep('combination', names(mcols(bins)))] <- paste0('combination.', unique(model$info$condition))
    }
    conditions <- sub('combination.', '', grep('combination', names(mcols(bins)), value=TRUE))
    if (is.null(combinations)) {
        comb.levels <- levels(mcols(bins)[,paste0('combination.', conditions[1])])
    } else {
        comb.levels <- combinations
    }
    if (is.null(marks)) {
        mark.levels <- unique(model$info$mark)
    } else {
        mark.levels <- marks
    }

    if (what %in% c('peaks','counts')) {
        ### Get fold enrichment
        enrich <- enrichmentAtAnnotation(model$bins, model$info, annotation, bp.around.annotation=bp.around.annotation, region=region, what=what, num.intervals=num.intervals, statistic=statistic)
    }
    ggplts <- list()
    maxfolds <- list()
    for (condition in conditions) {
        if (what == 'combinations') {
            ### Get fold enrichment
            bins$combination <- mcols(bins)[,paste0('combination.', condition)]
            enrich.cond <- enrichmentAtAnnotation(bins, model$info, annotation, bp.around.annotation=bp.around.annotation, region=region, what=what, num.intervals=num.intervals, statistic=statistic)
        } else {
            enrich.cond <- enrich
        }
        ### Prepare for plotting
        df <- reshape2::melt(enrich.cond)
        df$L1 <- factor(df$L1, levels=c('start','inside','end'))
        df <- rbind(df[df$L1 == 'start',], df[df$L1 == 'inside',], df[df$L1 == 'end',])
        if (length(region)>=2 & 'inside' %in% region) {
            df <- df[!(df$L1 == 'start' & df$lag > 0),]
            df <- df[!(df$L1 == 'end' & df$lag < 0),]
            df$position <- apply(data.frame(df$interval, df$lag), 1, max, na.rm = TRUE)
        } else if (length(region)==2 & ! 'inside' %in% region) {
            df <- df[!(df$L1 == 'start' & df$lag > 0),]
            df <- df[!(df$L1 == 'end' & df$lag < 0),]
            df$position <- df$lag
        } else if (length(region)==1) {
            df <- df[df$L1 == region,]
            df$position <- df$lag
        }
        if ('inside' %in% region) {
            df$position[df$L1 == 'end'] <- df$position[df$L1 == 'end'] + bp.around.annotation
        }
        df$position[df$L1 == 'inside'] <- df$position[df$L1 == 'inside'] * bp.around.annotation
        if (what == 'combinations') {
            df <- df[df$combination %in% comb.levels,]
        } else if (what %in% c('peaks','counts')) {
            df$mark <- sub("-.*", "", df$track)
            df <- df[df$mark %in% mark.levels, ]
            df$condition <- sapply(strsplit(as.character(df$track), '-'), '[[', 2)
            if (condition != "") {
                df <- df[df$condition == condition, ]
            }
        }

        ### Plot
        if (what == 'combinations') {
            ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='combination'), size=2)
            if (statistic == 'fold') {
                ggplt <- ggplt + ylab('fold enrichment')
                ggplt <- ggplt + geom_hline(yintercept=1, lty=2)
            } else if (statistic == 'fraction') {
                ggplt <- ggplt + ylab('fraction')
            }
            ggplt <- ggplt + scale_color_manual(values = getDistinctColors(length(unique(df$combination))))
        } else if (what == 'peaks') {
            ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='mark'), size=2)
            if (statistic == 'fold') {
                ggplt <- ggplt + ylab('fold enrichment')
                ggplt <- ggplt + geom_hline(yintercept=1, lty=2)
            } else if (statistic == 'fraction') {
                ggplt <- ggplt + ylab('fraction')
            }
            ggplt <- ggplt + scale_color_manual(values = getDistinctColors(length(unique(df$mark))))
        } else if (what == 'counts') {
            ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='track'), size=2) + ylab('RPKM')
            ggplt <- ggplt + scale_color_manual(values = getDistinctColors(length(unique(df$track))))
        }
        ggplt <- ggplt + theme_bw() + xlab('distance from annotation in [bp]')
        if (length(region)>=2 & 'inside' %in% region) {
            breaks <- c(c(-1, -0.5, 0, 0.5, 1, 1.5, 2) * bp.around.annotation)
            labels <- c(-bp.around.annotation, -bp.around.annotation/2, '0%', '50%', '100%', bp.around.annotation/2, bp.around.annotation)
            ggplt <- ggplt + scale_x_continuous(breaks=breaks, labels=labels)
        }
        maxfolds[[condition]] <- max(df$value, na.rm=TRUE)
        ggplts[[condition]] <- ggplt
    }
    maxfolds <- unlist(maxfolds)
    if (statistic == 'fraction' & what %in% c('combinations','peaks')) {
        ggplts <- lapply(ggplts, function(ggplt) { ggplt + scale_y_continuous(limits=c(0,1)) })
    } else {
        ggplts <- lapply(ggplts, function(ggplt) { ggplt + scale_y_continuous(limits=c(0,max(maxfolds, na.rm=TRUE)*1.1)) })
    }
    if (class(model) == class.multivariate.model) {
        return(ggplts[[1]])
    } else if (class(model) == class.combined.multivariate.model) {
        return(ggplts)
    }
    
}


#' Enrichment of (combinatorial) states for genomic annotations
#'
#' The function calculates the enrichment of a genomic feature with peaks or combinatorial states. Input is a \code{\link{multimodel}} object (containing the peak calls and combinatorial states) and a \code{\link{GRanges}} object containing the annotation of interest (e.g. transcription start sites or genes).
#'
#' @author Aaron Taudt
#' @param bins The \code{$bins} entry from a \code{\link{multimodel}} or \code{\link{combinedMultimodel}} object.
#' @param info The \code{$info} entry from a \code{\link{multimodel}} or \code{\link{combinedMultimodel}} object.
#' @param annotation A \code{\link{GRanges}} object with the annotation of interest.
#' @param bp.around.annotation An integer specifying the number of basepairs up- and downstream of the annotation for which the enrichment will be calculated.
#' @param region A combination of \code{c('start','inside','end')} specifying the region of the annotation for which the enrichment will be calculated. Select \code{'start'} if you have a point-sized annotation like transcription start sites. Select \code{c('start','inside','end')} if you have long annotations like genes.
#' @param what One of \code{c('combinations','peaks','counts','transitions')} specifying which statistic to calculate.
#' @param num.intervals Number of intervals for enrichment 'inside' of annotation.
#' @param statistic The statistic to calculate. Either 'fold' for fold enrichments or 'fraction' for fraction of bins falling into the annotation.
#' @return A \code{list()} containing \code{data.frame()}s for enrichment of combinatorial states and binary states at the start, end and inside of the annotation.
#' @importFrom S4Vectors as.factor subjectHits queryHits
enrichmentAtAnnotation <- function(bins, info, annotation, bp.around.annotation=10000, region=c('start','inside','end'), what=c('combinations','peaks','counts'), num.intervals=21, statistic='fold') {

    ## Check user input
    if ((!what %in% c('combinations','peaks','counts')) | length(what) > 1) {
        stop("argument 'what' must be one of c('combinations','peaks','counts')")
    }
  
    ## Variables
    binsize <- width(bins)[1]
    lag <- round(bp.around.annotation/binsize)
    enrich <- list()
    enrich$combinations <- list()
    enrich$peaks <- list()
    enrich$counts <- list()
    info.dedup <- info[!duplicated(paste0(info$mark, info$condition)), ]

    ## Get combinatorial and binary states
    combinations <- bins$combination
    tcombinations <- table(combinations)
    if ('peaks' %in% what) {
        binstates <- dec2bin(bins$state, colnames=info$ID)
        # Remove replicates
        binstates <- binstates[ ,info.dedup$ID]
        colsums.binstates <- colSums(binstates)
    }
    if ('counts' %in% what) {
        counts <- bins$counts
    }

    ### Calculating enrichment inside of annotation ###
    if ('inside' %in% region) {
        ptm <- startTimedMessage("Enrichment inside of annotations ...")

        intervals <- seq(from=0, to=1, length.out=num.intervals+1)
        widths.annotation <- width(annotation) - 1
        annotation.1bp <- resize(annotation, 1, fix='start')
        # Initialize arrays
        if ('peaks' %in% what) binstates.inside <- array(dim=c(num.intervals+1, length(info.dedup$ID)), dimnames=list(interval=intervals, track=info.dedup$ID))
        if ('combinations' %in% what) combinations.inside <- array(dim=c(num.intervals+1, length(levels(bins$combination))), dimnames=list(interval=intervals, combination=levels(bins$combination)))
        if ('counts' %in% what) counts.inside <- array(dim=c(num.intervals+1, length(info$ID)), dimnames=list(interval=intervals, track=info$ID))

        for (interval in intervals) {
            shift <- widths.annotation * interval * c(1,-1,1)[as.integer(strand(annotation))]
            shifted.starts <- start(annotation.1bp) + shift
            annotation.shifted <- GRanges(seqnames = seqnames(annotation.1bp), ranges = IRanges(start = shifted.starts, end = shifted.starts), strand = strand(annotation.1bp))
            # Get bins that overlap the shifted annotation
            index.inside.plus <- findOverlaps(annotation.shifted[strand(annotation.shifted)=='+' | strand(annotation.shifted)=='*'], bins, select="first")
            index.inside.minus <- findOverlaps(annotation.shifted[strand(annotation.shifted)=='-'], bins, select="last")
            index.inside.plus <- index.inside.plus[!is.na(index.inside.plus)]
            index.inside.minus <- index.inside.minus[!is.na(index.inside.minus)]
            index <- c(index.inside.plus, index.inside.minus)
            index <- index[index>0 & index<=length(bins)] # index could cross chromosome boundaries, but we risk it
            if ('peaks' %in% what) {
                if (statistic == 'fraction') {
                    binstates.inside[as.character(interval),] <- colSums(binstates[index,]) / length(index) # or colMeans
                } else if (statistic == 'fold') {
                    binstates.inside[as.character(interval),] <- colSums(binstates[index,]) / length(index) / colsums.binstates * length(bins)
                }
            }
            if ('combinations' %in% what) {
                if (statistic == 'fraction') {
                    fold <- table(combinations[index]) / length(index)
                } else if (statistic == 'fold') {
                    fold <- table(combinations[index]) / length(index) / tcombinations * length(bins) # fold enrichment
                }
                fold[is.na(fold)] <- 0
                combinations.inside[as.character(interval),] <- fold
            }
            if ('counts' %in% what) {
                counts.inside[as.character(interval),] <- colMeans(counts[index,])
            }
        }
        if ('peaks' %in% what) {
            enrich$peaks$inside <- binstates.inside
        }
        if ('combinations' %in% what) {
            enrich$combinations$inside <- combinations.inside
        }
        if ('counts' %in% what) {
            enrich$counts$inside <- counts.inside
        }
        stopTimedMessage(ptm)
    }

    ### 10000 bp before annotation ###
    if ('start' %in% region) {
        ptm <- startTimedMessage("Enrichment ",bp.around.annotation,"bp before annotations")
        # Get bins that overlap the start of annotation
        index.start.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], bins, select="first")
        index.start.minus <- findOverlaps(annotation[strand(annotation)=='-'], bins, select="last")
        index.start.plus <- index.start.plus[!is.na(index.start.plus)]
        index.start.minus <- index.start.minus[!is.na(index.start.minus)]
        # Occurrences at every bin position relative to feature
        if ('peaks' %in% what) binstates.start <- array(dim=c(length(-lag:lag), length(info.dedup$ID)), dimnames=list(lag=-lag:lag, track=info.dedup$ID))
        if ('combinations' %in% what) combinations.start <- array(dim=c(length(-lag:lag), length(levels(bins$combination))), dimnames=list(lag=-lag:lag, combination=levels(bins$combination)))
        if ('counts' %in% what) counts.start <- array(dim=c(length(-lag:lag), length(info$ID)), dimnames=list(lag=-lag:lag, track=info$ID))
        for (ilag in -lag:lag) {
            index <- c(index.start.plus+ilag, index.start.minus-ilag)
            index <- index[index>0 & index<=length(bins)]
            if ('peaks' %in% what) {
                if (statistic == 'fraction') {
                    binstates.start[as.character(ilag),] <- colSums(binstates[index,]) / length(index)
                } else if (statistic == 'fold') {
                    binstates.start[as.character(ilag),] <- colSums(binstates[index,]) / length(index) / colsums.binstates * length(bins)
                }
            }
            if ('combinations' %in% what) {
                if (statistic == 'fraction') {
                    fold <- table(combinations[index]) / length(index)
                } else if (statistic == 'fold') {
                    fold <- table(combinations[index]) / length(index) / tcombinations * length(bins) # fold enrichment
                }
                fold[is.na(fold)] <- 0
                combinations.start[as.character(ilag),] <- fold
            }
            if ('counts' %in% what) {
                counts.start[as.character(ilag),] <- colMeans(counts[index,])
            }
        }
        if ('peaks' %in% what) {
            rownames(binstates.start) <- as.numeric(rownames(binstates.start)) * binsize
            enrich$peaks$start <- binstates.start
        }
        if ('combinations' %in% what) {
            rownames(combinations.start) <- as.numeric(rownames(combinations.start)) * binsize
            enrich$combinations$start <- combinations.start
        }
        if ('counts' %in% what) {
            rownames(counts.start) <- as.numeric(rownames(counts.start)) * binsize
            enrich$counts$start <- counts.start
        }
        stopTimedMessage(ptm)
    }

    ### 10000 bp after annotation ###
    if ('end' %in% region) {
        ptm <- startTimedMessage("Enrichment ",bp.around.annotation,"bp after annotations")
        # Get bins that overlap the end of annotation
        index.end.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], bins, select="last")
        index.end.minus <- findOverlaps(annotation[strand(annotation)=='-'], bins, select="first")
        index.end.plus <- index.end.plus[!is.na(index.end.plus)]
        index.end.minus <- index.end.minus[!is.na(index.end.minus)]
        # Occurrences at every bin position relative to feature
        if ('peaks' %in% what) binstates.end <- array(dim=c(length(-lag:lag), length(info.dedup$ID)), dimnames=list(lag=-lag:lag, track=info.dedup$ID))
        if ('combinations' %in% what) combinations.end <- array(dim=c(length(-lag:lag), length(levels(bins$combination))), dimnames=list(lag=-lag:lag, combination=levels(bins$combination)))
        if ('counts' %in% what) counts.end <- array(dim=c(length(-lag:lag), length(info$ID)), dimnames=list(lag=-lag:lag, track=info$ID))
        for (ilag in -lag:lag) {
            index <- c(index.end.plus+ilag, index.end.minus-ilag)
            index <- index[index>0 & index<=length(bins)]
            if ('peaks' %in% what) {
                if (statistic == 'fraction') {
                    binstates.end[as.character(ilag),] <- colSums(binstates[index,]) / length(index)
                } else if (statistic == 'fold') {
                    binstates.end[as.character(ilag),] <- colSums(binstates[index,]) / length(index) / colsums.binstates * length(bins)
                }
            }
            if ('combinations' %in% what) {
                if (statistic == 'fraction') {
                    fold <- table(combinations[index]) / length(index)
                } else if (statistic == 'fold') {
                    fold <- table(combinations[index]) / length(index) / tcombinations * length(bins) # fold enrichment
                }
                fold[is.na(fold)] <- 0
                combinations.end[as.character(ilag),] <- fold
            }
            if ('counts' %in% what) {
                counts.end[as.character(ilag),] <- colMeans(counts[index,])
            }
        }
        if ('peaks' %in% what) {
            rownames(binstates.end) <- as.numeric(rownames(binstates.end)) * binsize
            enrich$peaks$end <- binstates.end
        }
        if ('combinations' %in% what) {
            rownames(combinations.end) <- as.numeric(rownames(combinations.end)) * binsize
            enrich$combinations$end <- combinations.end
        }
        if ('counts' %in% what) {
            rownames(counts.end) <- as.numeric(rownames(counts.end)) * binsize
            enrich$counts$end <- counts.end
        }
        stopTimedMessage(ptm)
    }

    return(enrich[[what]])

}



