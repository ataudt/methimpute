#' Get state colors
#'
#' Get the colors that are used for plotting.
#'
#' @param states Any combination of \code{c("unmethylated","methylated","heterozygous", "total")}.
#' @return A character vector with colors.
#' @seealso \code{\link{plotting}}
#' @export
#'@examples
#'cols <- getStateColors()
#'pie(1:length(cols), col=cols, labels=names(cols))
getStateColors <- function(states=NULL) {
    state.colors <- c("Background" = "gray", "UNmethylated" = "red","Methylated" = "blue","Heterozygous" = "green", "total" = "black", "background" = "gray", "signal" = "red")
    if (is.null(states)) {
        return(state.colors)
    } else {
        states <- as.character(states)
        states.other <- setdiff(states, names(state.colors))
        if (length(states.other) > 0) {
            state.colors.other <- getDistinctColors(states.other)
            state.colors <- c(state.colors, state.colors.other)
        }
        return(state.colors[states])
    }
}
 

#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges}} object.
#' @return The input \code{\link{GRanges}} with two additional metadata columns 'start.genome' and 'end.genome'.
transCoord <- function(gr) {
  cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
  cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- seqlevels(gr)
  gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  return(gr)
}


#' Plot a count histogram
#' 
#' Plot a histogram of count values and fitted distributions.
#' 
#' @export
plotHistogram <- function(model, total.counts, binwidth=1) {
  
    ## Get cross section at total counts ##
    counts <- model$data$counts
    mask.crosssec <- counts[,'total'] == total.counts
    counts.methylated <- counts[mask.crosssec,'methylated']
    
    contexts <- intersect(levels(model$data$context), unique(model$data$context))
    ggplts <- list()
    for (context in contexts) {
        mask.context <- model$data$context[mask.crosssec] == context
        ## Histogram ##
        ggplt <- ggplot(data.frame(counts=counts.methylated[mask.context])) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=binwidth, color='black', fill='white')
        ggplt <- ggplt + theme_bw()
        ggplt <- ggplt + xlab(paste0("methylated counts\nat positions with total counts = ", total.counts))
        ggplt <- ggplt + coord_cartesian(xlim=c(0,total.counts))
        
        ## Fits ##
        if (!is.null(model$params$emissionParams)) {
            x <- seq(0, total.counts, by = binwidth)
            distr <- list(x=x)
            for (i1 in 1:length(model$params$emissionParams)) {
                e <- model$params$emissionParams[[i1]]
                distr[[names(model$params$emissionParams)[i1]]] <- model$params$weights[[context]][i1] * dbinom(x, size=total.counts, prob=e[context,'prob'])
            }
            distr <- as.data.frame(distr)
            distr <- reshape2::melt(distr, id.vars='x', variable.name='components')
            distr$components <- sub('^X', '', distr$components)
            distr$components <- factor(distr$components, levels=levels(model$data$state))
            ggplt <- ggplt + geom_line(data=distr, mapping=aes_string(x='x', y='value', col='components'))
            
            ## Make legend
            lprobs <- round(sapply(model$params$emissionParams, function(x) { x[context,'prob'] }), 4)
            lweights <- round(model$params$weights[[context]], 2)
            legend <- paste0(names(model$params$emissionParams), ", prob=", lprobs, ", weight=", lweights)
            ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(rownames(model$params$emissionParams)), labels=legend) + theme(legend.position=c(0.5,0.99), legend.justification=c(0.5,1))
        }
        ggplt <- ggplt + ggtitle(paste0("Context ", context))
        ggplts[[context]] <- ggplt
    }
    cowplt <- cowplot::plot_grid(plotlist = ggplts, align='v', ncol=1)
        
    return(cowplt)
      
}


#' Plot a ratio histogram
#'
#' Plot a histogram of ratio values and fitted distributions.
#'
plotRatioHistogram <- function(model, num.intervals=20) {

    ## Plot histogram
    ggplt <- ggplot(data.frame(ratio=model$data$ratio)) + geom_histogram(aes_string(x='ratio', y='..density..'), breaks=seq(0, 1, length.out = num.intervals+1), color='black', fill='white') + coord_cartesian(xlim=c(0,1)) + xlab("ratio")

    ## Add distributions
    x <- c(seq(0, 0.1, by=0.001), seq(0.1, 0.9, by=0.01), seq(0.9, 1, by=0.001))
    distr <- list(x=x)
    p <- model$params
    for (irow in 1:nrow(p$emissionParams)) {
        e <- p$emissionParams
        distr[[rownames(p$emissionParams)[irow]]] <- p$weights[irow] * dbeta(x, e[irow,'a'], e[irow,'b'])
    }
    distr <- as.data.frame(distr)
    distr$total <- rowSums(distr[,2:4])
    distr <- reshape2::melt(distr, id.vars='x', variable.name='components')
    ggplt <- ggplt + geom_line(data=distr, mapping=aes_string(x='x', y='value', col='components'))
    
    ## Make legend
    lweights <- round(p$weights, 2)
    legend <- paste0(rownames(p$emissionParams), ", weight=", lweights)
    legend <- c(legend, paste0('Total, weight=1'))
    ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(c(rownames(p$emissionParams), 'total')), labels=legend) + theme(legend.position=c(0.5,1), legend.justification=c(0.5,1))
    return(ggplt)
}


plotRatioBoxplot <- function(model) {
  
    df <- data.frame(state=model$data$state, ratio=model$data$ratio)
    ggplt <- ggplot(df) + geom_boxplot(aes_string(x='state', y='ratio'))
    return(ggplt)
    
}


#' Plot a count histogram
#' 
#' Plot a histogram of count values and fitted distributions.
#' 
plotHistogram2 <- function(model, binwidth=10) {
    
    ## Assign variables
    counts <- model$data$observable
    maxcounts <- max(counts)
    
    ## Plot histogram
    ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=binwidth, color='black', fill='white') + xlab("counts")
    ggplt <- ggplt + coord_cartesian(xlim=c(0,quantile(counts, 0.995)))
    ggplt <- ggplt + theme_bw()
    
    ## Add distributions
    if (!is.null(model$params$emissionParams)) {
        x <- seq(0, maxcounts, by = binwidth)
        distr <- list(x=x)
        for (irow in 1:nrow(model$params$emissionParams)) {
            e <- model$params$emissionParams
            if (e$type[irow] == 'delta') {
                distr[[rownames(model$params$emissionParams)[irow]]] <- rep(0, length(x))
                distr[[rownames(model$params$emissionParams)[irow]]][1] <- model$params$weights[irow]
            } else if (e$type[irow] == 'dnbinom') {
                distr[[rownames(model$params$emissionParams)[irow]]] <- model$params$weights[irow] * dnbinom(x, size=e[irow,'size'], prob=e[irow,'prob'])
            }
        }
        distr <- as.data.frame(distr)
        distr$total <- rowSums(distr[,2:(1+nrow(model$params$emissionParams))])
        distr <- reshape2::melt(distr, id.vars='x', variable.name='components')
        distr$components <- sub('^X', '', distr$components)
        distr$components <- factor(distr$components, levels=c(levels(model$data$state), 'total'))
        ggplt <- ggplt + geom_line(data=distr, mapping=aes_string(x='x', y='value', col='components'))
        
        ## Make legend
        lmeans <- round(model$params$emissionParams[,'mu'], 2)
        lvars <- round(model$params$emissionParams[,'var'], 2)
        lweights <- round(model$params$weights, 2)
        legend <- paste0(rownames(model$params$emissionParams), ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
        legend <- c(legend, paste0('total, mean=', round(mean(counts),2), ', var=', round(var(counts),2)))
        ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(c(rownames(model$params$emissionParams),'total')), labels=legend) + theme(legend.position=c(1,1), legend.justification=c(1,1))
    }
    
    return(ggplt)
}


plotBoxplot <- function(model, datapoints=1000) {
  
    df <- data.frame(state=model$data$state, observable=model$data$observable)
    if (datapoints < nrow(df)) {
        df <- df[sample(1:nrow(df), datapoints, replace = FALSE), ]
    }
    df <- suppressMessages( reshape2::melt(df) )
    names(df) <- c('state', 'status', 'observable')
    ggplt <- ggplot(df) + geom_boxplot(aes_string(x='state', y='observable', fill='status'))
    return(ggplt)
    
}


plotBoxplotRatio <- function(model) {
  
    df <- data.frame(state=model$data$state, ratio=model$data$ratio)
    ggplt <- ggplot(df) + geom_boxplot(aes_string(x='state', y='ratio'))
    return(ggplt)
    
}


plotScatter <- function(model, datapoints=1000) {
  
    contexts <- intersect(levels(model$data$context), unique(model$data$context))
    ggplts <- list()
    limits <- list()
    for (context in contexts) {
        data <- model$data[model$data$context == context]
        ## Find sensible limits
        xmax <- quantile(data$counts[,'total']-data$counts[,'methylated'], 0.99)
        ymax <- quantile(data$counts[,'methylated'], 0.99)
        limits[[context]] <- c(xmax, ymax)
        df <- data.frame(state=data$state, unmethylated=data$counts[,'total']-data$counts[,'methylated'], methylated=data$counts[,'methylated'])
        if (datapoints < nrow(df)) {
            df <- df[sample(1:nrow(df), datapoints, replace = FALSE), ]
        }
        
        ggplt <- ggplot(df, aes_string(x='methylated', y='unmethylated', col='state'))
        ggplt <- ggplt + geom_point(alpha=0.3)
        ggplt <- ggplt + coord_cartesian(xlim=c(0,xmax), ylim=c(0,ymax))
        ggplt <- ggplt + theme_bw()
        
        ## Legend
        lweights <- round(model$params$weights[[context]], 2)
        legend <- paste0(names(model$params$weights[[context]]), ", weight=", lweights)
        ggplt <- ggplt + scale_color_manual(values=getStateColors(names(model$params$weights[[context]])), labels=legend)
        ggplt <- ggplt + ggtitle(paste0("Context ", context))
        ggplt <- ggplt + theme(legend.position=c(1,1), legend.justification=c(1,1))
        ggplts[[context]] <- ggplt
    }
    xmax <- max(sapply(limits, '[[', 1))
    ymax <- max(sapply(limits, '[[', 2))
    ggplts <- lapply(ggplts, function(ggplt) { ggplt <- ggplt + coord_cartesian(xlim=c(0,xmax), ylim=c(0,ymax)); ggplt })
    cowplt <- cowplot::plot_grid(plotlist = ggplts, align = 'h', nrow=1)
    
    return(cowplt)
    
}


# =================================================================
# Plot a clustered heatmap of state calls
# =================================================================
#' Genome wide heatmap of CNV-state
#'
#' Plot a genome wide heatmap of copy number variation state. This heatmap is best plotted to file, because in most cases it will be too big for cleanly plotting it to screen.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the same length as \code{hmms}. If \code{NULL} the IDs from the \code{\link{aneuHMM}} objects will be used.
#' @param classes A character vector with the classification of the elements on the y-axis. The vector must have the same length as \code{hmms}.
#' @param reorder.by.class If \code{TRUE}, the dendrogram will be reordered to display similar classes next to each other.
#' @param classes.color A (named) vector with colors that are used to distinguish \code{classes}. Names must correspond to the unique elements in \code{classes}.
#' @param file A PDF file to which the heatmap will be plotted.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @param range A \code{\link[GenomicRanges]{GRanges}} object to subset the range for which the heatmap will be plotted.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
#' @importFrom stats as.dendrogram
#' @importFrom ggdendro dendro_data theme_dendro
#' @export
#'@examples
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
#'## Plot a clustered heatmap
#'classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files)))
#'labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files)))
#'heatmapGenomewide(c(lung.files, liver.files), ylabels=labels, classes=classes,
#'                  classes.color=c('blue','red'))
#'
heatmapGenomewide <- function(models, ylabels=NULL, classes=NULL, reorder.by.class=TRUE, classes.color=NULL, file=NULL, cluster=TRUE, range=NULL) {

  	## Check user input
  	if (!is.null(ylabels)) {
    		if (length(ylabels) != length(models)) {
      			stop("length(ylabels) must equal length(models)")
    		}
  	}
  	if (!is.null(classes)) {
    		if (length(classes) != length(models)) {
      			stop("length(classes) must equal length(models)")
    		}
  	}
  	if (length(classes.color)!=length(unique(classes))) {
    		stop("'classes.color' must have the same length as unique(classes)")
  	}
  	if (is.null(names(classes.color))) {
    		names(classes.color) <- unique(classes)
  	}
  	if (!setequal(names(classes.color), unique(classes))) {
    		stop("The names of 'classes.color' must be equal to the unique elements in 'classes'")
  	}

  	## Load the files
  	models <- loadFromFiles(models, check.class='binnedMethylome')
  
  	## Dataframe with IDs, ylabels and classes
  	class.data <- data.frame(ID=sapply(models,'[[','ID'))
  	class.data$ID <- factor(class.data$ID, levels=class.data$ID)
  	if (is.null(ylabels)) {
    	  class.data$ylabel <- as.character(class.data$ID)
  	} else {
      	class.data$ylabel <- as.character(ylabels)
  	}
  	class.data$class <- classes
  	
  	## Mapping to match ID with ylabel
  	mapping <- class.data$ylabel
  	names(mapping) <- class.data$ID

  	## Get segments and SCE coordinates
  	if (reorder.by.class) {
      	temp <- getSegments(models, cluster=cluster, classes=classes)
  	} else {
      	temp <- getSegments(models, cluster=cluster)
  	}
  	segments.list <- temp$segments
  	hc <- temp$clustering
  	if (cluster) {
    		models <- models[hc$order]
    		class.data <- class.data[hc$order,]
    		class.data$ID <- factor(class.data$ID, levels=class.data$ID)
  	}
  	
  	## Subsetting coordinates
  	if (!is.null(range)) {
  	    range <- transCoord(range)
  	    segments.list <- endoapply(segments.list, function(x) { x <- subsetByOverlaps(x, range); x <- keepSeqlevels(x, intersect(seqlevels(x), unique(seqnames(x)))); x })
  	}

  	## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
  	ptm <- startTimedMessage("Transforming coordinates ...")
  	segments.list <- endoapply(segments.list, transCoord)
  	stopTimedMessage(ptm)

  	## Data.frame for plotting
  	ptm <- startTimedMessage("Making the plot ...")
  	df <- list()
  	for (i1 in 1:length(segments.list)) {
    		df[[length(df)+1]] <- data.frame(start=segments.list[[i1]]$start.genome, end=segments.list[[i1]]$end.genome, seqnames=seqnames(segments.list[[i1]]), ID=names(segments.list)[i1], state=segments.list[[i1]]$state)
  	}
  	df <- do.call(rbind, df)
  	df$ID <- factor(df$ID, levels=levels(class.data$ID))
  	df$ylabel <- mapping[as.character(df$ID)]
	
  	# Chromosome lines
  	gr.chromlines <- GRanges(seqnames=seqlevels(segments.list[[1]]), ranges=IRanges(start=seqlengths(segments.list[[1]]), end=seqlengths(segments.list[[1]])))
  	seqlengths(gr.chromlines) <- seqlengths(segments.list[[1]])
  	gr.chromlines <- transCoord(gr.chromlines)
  	if (!is.null(range)) {
  	    gr.chromlines <- subsetByOverlaps(gr.chromlines, range)
  	}
  	if (length(gr.chromlines) > 0) {
      	label.pos <- round( c(0, gr.chromlines$start.genome[1:(length(gr.chromlines)-1)]) + 0.5 * seqlengths(gr.chromlines) )
      	names(label.pos) <- seqlevels(gr.chromlines)
      	df.chroms <- data.frame(y=c(0, gr.chromlines$start.genome), x=1, xend=length(segments.list))
  	}

  	### Plot ###
  	pltlist <- list()
  	widths <- vector()

  	## Prepare the plot
  	df$x <- as.numeric(df$ID) # transform all x-coordiantes to numeric because factors and numerics get selected different margins
  	ggplt <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='x', col='state'), size=5)
  	ggplt <- ggplt + scale_x_continuous(name="sample", breaks=1:length(unique(df$ylabel)), labels=unique(df$ylabel))
  	ggplt <- ggplt + scale_color_manual(values=getStateColors(levels(df$state)))
  	ggplt <- ggplt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20), axis.line=element_blank(), axis.title.x=element_blank())
  	# ggplt <- ggplt + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
  	if (length(gr.chromlines) > 0) {
      	ggplt <- ggplt + scale_y_continuous(breaks=label.pos, labels=names(label.pos))
      	ggplt <- ggplt + geom_segment(aes_string(x='x', xend='xend', y='y', yend='y'), data=df.chroms, col='black')
  	}
  	ggplt <- ggplt + coord_flip()
  	# width.heatmap <- sum(as.numeric(seqlengths(hmms[[1]]$data))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
  	width.heatmap <- 150
  	height <- length(models) * 0.5
  	pltlist[['heatmap']] <- ggplt
  	widths['heatmap'] <- width.heatmap
  	## Make the classification bar
  	if (!is.null(classes)) {
    		width.classes <- 5
      	class.data$x <- as.numeric(class.data$ID)  # transform all x-coordiantes to numeric because factors and numerics get selected different margins
    		ggclass <- ggplot(class.data) + geom_linerange(aes_string(ymin=0, ymax=1, x='x', col='class'), size=5) + guides(col=FALSE) + xlab("class")
      	ggclass <- ggclass + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title.x=element_blank())
    		ggclass <- ggclass + coord_flip()
    		if (!is.null(classes.color)) {
      			ggclass <- ggclass + scale_color_manual(breaks=names(classes.color), values=classes.color)
    		}
    		pltlist[['classbar']] <- ggclass
    		widths['classbar'] <- width.classes
  	}
  	## Prepare the dendrogram
  	if (!is.null(hc)) {
    		dhc <- stats::as.dendrogram(hc)
    		ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    		ggdndr <- ggplot(ddata$segments) + geom_segment(aes_string(x='x', xend='xend', y='y', yend='yend')) + scale_y_reverse()
    		ggdndr <- ggdndr + coord_flip()
    		ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
    		width.dendro <- 20
    		pltlist[['dendro']] <- ggdndr
    		widths['dendro'] <- width.dendro
  	}
  	cowplt <- cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths=rev(widths))
  	stopTimedMessage(ptm)

  	## Plot to file
  	if (!is.null(file)) {
    		ptm <- startTimedMessage("Plotting to file ",file," ...")
    		ggsave(file, cowplt, width=sum(widths)/2.54, height=height/2.54, limitsize=FALSE)
    		stopTimedMessage(ptm)
  	} else {
    		return(cowplt)
  	}

}


#' Heatmap of transition probabilities
#'
#' Plot a heatmap of transition probabilities for a \code{\link{multiHMM}} model.
#'
#' @param model A \code{\link{multiHMM}} object or file that contains such an object.
#' @param order Set to \code{TRUE} if you want to order the heatmap.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @importFrom reshape2 melt
#' @seealso \code{\link{plotting}}
#' @export
heatmapTransitionProbs <- function(model, order=FALSE) {
  
    if (is.list(model$params$transProbs)) {
        As <- model$params$transProbs
    } else if (is.matrix(model$params$transProbs)) {
        As <- list(transitions=model$params$transProbs)
    }

    ggplts <- list()
    for (i1 in 1:length(As)) {
        # model <- suppressMessages( loadFromFiles(model, check.class=class.univariate.hmm)[[1]] )
        A <- reshape2::melt(As[[i1]], varnames=c('from','to'), value.name='prob')
        if (order) {
            stateorder <- stateorderByTransition(model$params$transProbs)
            A$from <- factor(A$from, levels=stateorder)
            A$to <- factor(A$to, levels=stateorder)
        } else {
            A$from <- factor(A$from)
            A$to <- factor(A$to)
        }
        ggplt <- ggplot(data=A) + geom_tile(aes_string(x='to', y='from', fill='prob')) + scale_fill_gradient(low="white", high="blue", limits=c(0,1)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
        ggplt <- ggplt + ggtitle(names(model$params$transProbs)[i1])
        ggplts[[names(model$params$transProbs)[i1]]] <- ggplt
    }
    if (length(ggplts) > 1) {
        # ggplts <- insertNULL(ggplts)
        ggplts <- cowplot::plot_grid(plotlist = ggplts, align='hv')
    }
    
    return(ggplts)
  
}

insertNULL <- function(plotlist) {
    n <- 0.5 * ( sqrt(8*length(plotlist)+1) - 1 )
    plotlist2 <- list()
    i3 <- 1
    for (i1 in 1:n) {
        if (i1 > 1) {
            for (i2 in 1:(n-i1)) {
                plotlist2[length(plotlist2)+1] <- list(NULL)
            }
        }
        for (i2 in i1:n) {
            plotlist2[[length(plotlist2)+1]] <- plotlist[[i3]]
            i3 <- i3 + 1
        }
    }
    return(plotlist2)
}

#' Histograms faceted by state
#' 
#' Make histograms of variables in \code{bins} faceted by state in \code{segmentation}.
#' 
#' @param segmentation A segmentation with metadata column 'state'.
#' @param bins A \code{\link{binnedMethylome}}.
#' @param plotcol A character vector specifying the meta-data columns for which to produce the faceted plot.
#' @param binwidth Bin width for the histogram.
#' @return A \code{\link[ggplot2]{ggplot}}.
#' @importFrom reshape2 melt
plotStateHistograms <- function(segmentation, bins, plotcol, binwidth) {
  
    states <- levels(segmentation$state)
    segmentation.split <- split(segmentation, segmentation$state)
    mind <- findOverlaps(bins, segmentation.split, select='first')
    bins$state <- factor(names(segmentation.split)[mind], levels=states)
    
    df <- data.frame(state=bins$state, mcols(bins)[,plotcol])
    names(df)[2:ncol(df)] <- plotcol
    dfm <- reshape2::melt(df, id.vars = 'state', measure.vars = plotcol, variable.name = 'plotcol')
    ggplt <- ggplot(dfm) + geom_freqpoly(aes_string(x='value', y='..density..', col='plotcol'), binwidth = binwidth)
    ggplt <- ggplt + theme_bw()
    ggplt <- ggplt + facet_wrap(~ state)
    
    return(ggplt)
}


#' Scatter plots faceted by state
#' 
#' Make scatter plots of two variables in \code{bins} faceted by state in \code{segmentation}.
#' 
#' @param segmentation A segmentation with metadata column 'state'.
#' @param bins A \code{\link{binnedMethylome}}.
#' @param x A character specifying the meta-data column that should go on the x-axis.
#' @param y A character specifying the meta-data column that should go on the y-axis.
#' @param col A character specifying the meta-data column that should be used for coloring.
#' @return A \code{\link[ggplot2]{ggplot}}.
plotStateScatter <- function(segmentation, bins, x, y, col=NULL, datapoints=1000) {
  
    states <- levels(segmentation$state)
    segmentation.split <- split(segmentation, segmentation$state)
    mind <- findOverlaps(bins, segmentation.split, select='first')
    bins$seg.state <- factor(names(segmentation.split)[mind], levels=states)
    
    df <- data.frame(state = bins$seg.state, x = mcols(bins)[,x], y = mcols(bins)[,y], color=mcols(bins)[,col])
    df <- df[sample(1:nrow(df), size=datapoints, replace=FALSE),]
    if (!is.null(col)) {
        ggplt <- ggplot(df) + geom_jitter(aes_string(x='x', y='y', col='color'), alpha=0.1)
    } else {
        ggplt <- ggplot(df) + geom_jitter(aes_string(x='x', y='y'), alpha=0.1)
    }
    ggplt <- ggplt + theme_bw()
    ggplt <- ggplt + xlab(x) + ylab(y)
    ggplt <- ggplt + facet_wrap(~ state)
    ggplt <- ggplt + scale_color_manual(values=getDistinctColors(levels(df$color)))
    
    return(ggplt)
}


#' Plot convergence info
#' 
#' Plot the convergence of the probability parameter in different sequence contexts.
#' 
#' @param model A \code{\link{BinomialHMM}} object or file that contains such an object.
#' @return A \code{\link[ggplot2]{ggplot}}.
plotConvergence <- function(model) {
    
    info <- model$convergenceInfo$prefit
    if (is.null(info)) {
        info <- model$convergenceInfo
    }
    # ## Plot loglikelihood
    # df <- data.frame(iteration=0:(length(info$logliks)-1), loglik=info$logliks)
    # ggplt <- ggplot(df) + geom_line(aes_string(x='iteration', y='loglik')) + theme_bw()
    # ggplt <- ggplt + scale_x_continuous(breaks=df$iteration, limits=c(2, length(info$logliks)-1))
    
    ## Plot parameters
    df <- reshape2::melt(info$parameterInfo, value.name='prob')
    ggplt <- ggplot(df) + geom_line(aes_string(x='iteration', y='prob', col='state')) + theme_bw()
    ggplt <- ggplt + theme(panel.grid.minor.x = element_blank())
    ggplt <- ggplt + facet_wrap(~ context)
    ggplt <- ggplt + scale_y_continuous(limits=c(0,1))
    ggplt <- ggplt + scale_color_manual(values=getStateColors(unique(df$state)))
    return(ggplt)
}


plotEnrichment <- function(data, annotation, windowsize=100, insidewindows=20, range=1000, category.column=NULL) {
  
    ## Variables
    categories <- 'none'
    if (!is.null(category.column)) {
        categories <- levels(mcols(data)[,category.column])
    }
  
    overlaps <- array(NA, dim=c(length(levels(data$context)), length(categories), range%/%windowsize, 2, 2), dimnames=list(context=levels(data$context), category=categories, distance=1:(range%/%windowsize), what=c('mean', 'weight'), where=c('upstream', 'downstream')))
    ## Upstream and downstream annotation
    annotation.up <- resize(x = annotation, width = 1, fix = 'start')
    annotation.down <- resize(x = annotation, width = 1, fix = 'end')
    for (context in levels(data$context)) {
        message("Upstream and downstream counts for context ", context)
        data.context <- data[data$context == context]
        for (category in categories) {
            message("  category ", category)
            if (!is.null(category.column)) {
                data.category <- data.context[mcols(data.context)[,category.column] == category]
            } else {
                data.category <- data.context
            }
            for (i1 in 1:(range %/% windowsize)) {
                anno.up <- suppressWarnings( promoters(x = annotation.up, upstream = i1*windowsize, downstream = 0) )
                anno.up <- suppressWarnings( resize(x = anno.up, width = windowsize, fix = 'start') )
                ind <- findOverlaps(anno.up, data.category)
                overlaps[context, category, as.character(i1), 'mean', 'upstream'] <- mean(data.category$posteriorMeth[ind@to])
                overlaps[context, category, as.character(i1), 'weight', 'upstream'] <- length(ind@to)
                # overlaps[context, category, as.character(-i1*windowsize), 'mean']] <- mean(data.category$posteriorMeth[ind@to])
                # overlaps[context, category, as.character(-i1*windowsize), 'weight']] <- length(ind@to)
            }
            for (i1 in 1:(range %/% windowsize)) {
                anno.down <- suppressWarnings( promoters(x = annotation.down, upstream = 0, downstream = i1*windowsize) )
                anno.down <- suppressWarnings( resize(x = anno.down, width = windowsize, fix = 'end') )
                ind <- findOverlaps(anno.down, data.category)
                overlaps[context, category, as.character(i1), 'mean', 'downstream'] <- mean(data.category$posteriorMeth[ind@to])
                overlaps[context, category, as.character(i1), 'weight', 'downstream'] <- length(ind@to)
                # overlaps[context, category, as.character((i1-1)*windowsize), 'mean']] <- mean(data.category$posteriorMeth[ind@to])
                # overlaps[context, category, as.character((i1-1)*windowsize), 'weight']] <- length(ind@to)
            }
        }
    }

    ## Inside annotation
    ioverlaps <- array(NA, dim=c(length(levels(data$context)), length(categories), insidewindows, 2, 1), dimnames=list(context=levels(data$context), category=categories, distance=1:insidewindows, what=c('mean', 'weight'), where=c('inside')))
    mask.plus <- as.logical(strand(annotation) == '+' | strand(annotation) == '*')
    mask.minus <- !mask.plus
    widths <- width(annotation)
    for (context in levels(data$context)) {
        message("Inside counts for context ", context)
        data.context <- data[data$context == context]
        for (category in categories) {
            message("  category ", category)
            if (!is.null(category.column)) {
                data.category <- data.context[mcols(data.context)[,category.column] == category]
            } else {
                data.category <- data.context
            }
            for (i1 in 1:insidewindows) {
                anno.in <- annotation
                end(anno.in[mask.plus]) <- start(annotation[mask.plus]) + round(widths[mask.plus] * i1 / insidewindows)
                start(anno.in[mask.plus]) <- start(annotation[mask.plus]) + round(widths[mask.plus] * (i1-1) / insidewindows)
                start(anno.in[mask.minus]) <- end(annotation[mask.minus]) - round(widths[mask.minus] * i1 / insidewindows)
                end(anno.in[mask.minus]) <- end(annotation[mask.minus]) - round(widths[mask.minus] * (i1-1) / insidewindows)
                ind <- findOverlaps(anno.in, data.category)
                ioverlaps[context, category, as.character(i1), 'mean', 'inside'] <- mean(data.category$posteriorMeth[ind@to])
                ioverlaps[context, category, as.character(i1), 'weight', 'inside'] <- length(ind@to)
                # overlaps[context, category, as.character((i1-1)*windowsize), 'mean']] <- mean(data.category$posteriorMeth[ind@to])
                # overlaps[context, category, as.character((i1-1)*windowsize), 'weight']] <- length(ind@to)
            }
        }
    }
    
    dfo <- reshape2::melt(overlaps)
    df <- dfo[dfo$what == 'mean',]
    names(df)[6] <- 'mean'
    df$distance <- as.numeric(df$distance)
    # Adjust distances for plot
    df$distance[df$where == 'upstream'] <- -df$distance[df$where == 'upstream'] * windowsize
    df$distance[df$where == 'downstream'] <- (df$distance[df$where == 'downstream']-1) * windowsize + range
    
    dfo <- reshape2::melt(ioverlaps)
    idf <- dfo[dfo$what == 'mean',]
    names(idf)[6] <- 'mean'
    idf$distance <- as.numeric(idf$distance)
    # Adjust distances for plot
    idf$distance <- (idf$distance-1) * range / insidewindows
    
    breaks <- c(c(-1, -0.5, 0, 0.5, 1, 1.5, 2) * range)
    labels <- c(-range, -range/2, '0%', '50%', '100%', range/2, range)
    
    weights <- apply(overlaps, c(1,2,4), sum)[,,'weight']
    iweights <- apply(ioverlaps, c(1,2,4), sum)[,,'weight']
    df <- rbind(df, idf)
    if (!is.null(category.column)) {
        df$weight <- diag((weights + iweights)[df$context, df$category])
    } else {
        df$weight <- (weights + iweights)[df$context]
    }
    df$weight <- log(df$weight+1)
    
    ggplt <- ggplot(df) + geom_line(aes_string(x='distance', y='mean', col='context', alpha='weight'))
    ggplt <- ggplt + scale_alpha_continuous(name='log(weight+1)', range=c(0.5,1))
    ggplt <- ggplt + scale_x_continuous(breaks=breaks, labels=labels)
    ggplt <- ggplt + xlab('Distance from annotation') + ylab('Posterior of methylation')
    if (!is.null(category.column)) {
        ggplt <- ggplt + facet_wrap(~ category)
    }
    ggplt
    
    return(ggplt)
}









