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
    state.colors <- c("UNmethylated" = "red","Methylated" = "blue","Hemimethylated" = "green", "total" = "black", "background" = "gray", "signal" = "red")
    if (is.null(states)) {
        return(state.colors)
    } else {
        return(state.colors[states])
    }
}
 

#' Plot a ratio histogram
#'
#' Plot a histogram of ratio values and fitted distributions.
#'
plotHistogramRatio <- function(hmm, num.intervals=20) {

    ## Plot histogram
    ggplt <- ggplot(data.frame(ratio=hmm$data$ratio)) + geom_histogram(aes_string(x='ratio', y='..density..'), breaks=seq(0, 1, length.out = num.intervals+1), color='black', fill='white') + coord_cartesian(xlim=c(0,1)) + xlab("ratio")

    ## Add distributions
    x <- c(seq(0, 0.1, by=0.001), seq(0.1, 0.9, by=0.01), seq(0.9, 1, by=0.001))
    distr <- list(x=x)
    p <- hmm$params
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

#' Plot a count histogram
#' 
#' Plot a histogram of count values and fitted distributions.
#' 
plotHistogram <- function(model) {
    
    ## Assign variables
    counts <- model$data$observable
    maxcounts <- max(counts)
    
    ## Plot histogram
    ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + xlab("counts")
    ggplt <- ggplt + theme_bw()
    
    ## Add distributions
    if (!is.null(model$params$emissionParams)) {
        x <- seq(0, maxcounts)
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