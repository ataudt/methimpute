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
    state.colors <- c("Background" = "gray", "UNmethylated" = "red","Methylated" = "blue","Hemimethylated" = "green", "total" = "black", "background" = "gray", "signal" = "red")
    if (is.null(states)) {
        return(state.colors)
    } else {
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
plotHistogram <- function(model, binwidth=1) {
    
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
            distr[[rownames(model$params$emissionParams)[irow]]] <- model$params$weights[irow] * dnbinom(x, size=e[irow,'size'], prob=e[irow,'prob'])
        }
        distr <- as.data.frame(distr)
        distr$total <- rowSums(distr[,2:(1+nrow(model$params$emissionParams))])
        distr <- reshape2::melt(distr, id.vars='x', variable.name='components')
        ggplt <- ggplt + geom_line(data=distr, mapping=aes_string(x='x', y='value', col='components'))
        
        ## Make legend
        lmeans <- round(model$params$emissionParams[,'mu'], 2)
        lvars <- round(model$params$emissionParams[,'var'], 2)
        lweights <- round(model$params$weights, 2)
        legend <- paste0(rownames(model$params$emissionParams), ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
        legend <- c(legend, paste0('total, mean=', round(mean(counts),2), ', var=', round(var(counts),2)))
        ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(c('background','signal','total')), labels=legend) + theme(legend.position=c(1,1), legend.justification=c(1,1))
    }
    
    return(ggplt)
}


plotBoxplot <- function(model) {
  
    df <- data.frame(state=model$data$state, observable=model$data$observable)
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


plotScatter <- function(model, datapoints='all') {
  
    ## Find sensible limits
    xmax <- quantile(model$data$counts.unmethylated, 0.99)
    ymax <- quantile(model$data$counts.methylated, 0.99)
    df <- data.frame(state=model$data$state, unmeth=model$data$counts.unmethylated, meth=model$data$counts.methylated)
    if (is.numeric(datapoints)) {
        df <- df[sample(1:nrow(df), datapoints, replace = FALSE), ]
    }
    
    ggplt <- ggplot(df, aes_string(x='unmeth', y='meth', col='state')) + theme_bw()
    ggplt <- ggplt + geom_point()
    # ggplt <- ggplt + geom_density2d()
    ggplt <- ggplt + scale_color_manual(values=getStateColors(names(model$params$weights)))
    ggplt <- ggplt + xlab('unmethylated counts') + ylab('methylated counts')
    ggplt <- ggplt + coord_cartesian(xlim=c(0,xmax), ylim=c(0,ymax))
    
    return(ggplt)
    
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

