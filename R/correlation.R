#' Distance correlation
#' 
#' Compute the distance correlation from a \code{\link{methimputeData}} object.
#' 
#' @param data A \code{\link{methimputeData}} object.
#' @param distances An integer vector specifying the distances for which the correlation will be calculated.
#' @param separate.contexts A logical indicating whether contexts are treated separately. If set to \code{TRUE}, correlations will only be calculated between cytosines of the same context.
#' @return A list() with an array containing the correlation values and the corresponding \code{\link[ggplot2]{ggplot}}.
#' 
#' @export
#' @examples 
#'## Get some toy data
#'file <- system.file("data","arabidopsis_toydata.RData",
#'                     package="methimpute")
#'data <- get(load(file))
#'distcor <- distanceCorrelation(data)
#'print(distcor$plot)
#'
distanceCorrelation <- function(data, distances=0:50, separate.contexts=FALSE) {
    
    ## Contexts
    contexts <- intersect(levels(data$context), unique(data$context))
    
    ## Add meth.lvl column
    data$meth.lvl <- data$counts[,'methylated'] / data$counts[,'total']
    
    ### Add distance and transition context to bins ###
    data$distance <- addDistance(data, separate.contexts=separate.contexts)
    data$transitionContext <- addTransitionContext(data, separate.contexts=separate.contexts)

    ## Loop through distances
    ptm <- startTimedMessage("Calculating correlations\n")
    # cors <- list()
    cor.array <- array(NA, dim=c(length(contexts), length(contexts), length(distances), 2), dimnames=list(context=contexts, context=contexts, distance=distances, what=c('correlation','weight')))
    for (i1 in distances) {
        if (i1 == distances[1]) {
            message("  for distance ", i1, appendLF=FALSE)
        } else {
            message(", ", i1, appendLF=FALSE)
        }
        ind <- which(data$distance == i1)
        meth.lvl <- data$meth.lvl[ind]
        meth.lvl.shift <- data$meth.lvl[ind-1]
        # ## All contexts
        # cor <- tryCatch(cor(meth.lvl, meth.lvl.shift, use='complete.obs'), error = function(e) NA, warning = function(w) NA)
        # cors[[as.character(i1)]] <- data.frame(correlation = cor, weight = length(ind))
        ## Context specific
        context.transition <- data$transitionContext[ind]
        cor.matrix <- matrix(NA, ncol=length(contexts), nrow=length(contexts), dimnames=list(context=contexts, context=contexts))
        weight.matrix <- matrix(NA, ncol=length(contexts), nrow=length(contexts), dimnames=list(context=contexts, context=contexts))
        for (c1 in 1:length(contexts)) {
            for (c2 in c1:length(contexts)) {
                context <- paste0(contexts[c1], '-', contexts[c2])
                context.rev <- paste0(contexts[c2], '-', contexts[c1])
                mask <- context.transition == context | context.transition == context.rev
                ratio1 <- meth.lvl[mask]
                ratio2 <- meth.lvl.shift[mask]
                cor.matrix[contexts[c1], contexts[c2]] <- tryCatch(cor(ratio1, ratio2, use='complete.obs'), error = function(e) NA, warning = function(w) NA)
                cor.matrix[contexts[c2], contexts[c1]] <- cor.matrix[contexts[c1], contexts[c2]]
                weight.matrix[contexts[c1], contexts[c2]] <- length(ratio1)
                weight.matrix[contexts[c2], contexts[c1]] <- weight.matrix[contexts[c1], contexts[c2]]
            }
        }
        cor.array[,,as.character(i1),'correlation'] <- cor.matrix
        cor.array[,,as.character(i1),'weight'] <- weight.matrix
    }
    message("\nFinished calculating correlations in", appendLF = FALSE)
    stopTimedMessage(ptm)
    
    ## Context correlation plots
    maxweights <- numeric()
    dfs <- list()
    for (c1 in 1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            context.transition <- paste0(contexts[c1], '-', contexts[c2])
            if (c1 <= c2) {
                df <- data.frame(distance = distances, correlation = cor.array[c1,c2,,'correlation'], weight = cor.array[c1,c2,,'weight'], from = contexts[c1], to = contexts[c2])
                df$logweight <- log(df$weight+1)
                maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
                dfs[[context.transition]] <- df
            }
        }
    }
    maxweight <- max(maxweights, na.rm = TRUE)
    miny <- min(cor.array, na.rm = TRUE)
    
    ## Plot correlation
    df <- do.call(rbind, dfs)
    ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
    ggplt <- ggplt + xlab('distance in [bp]')
    ggplt <- ggplt + facet_grid(from ~ to)
    if (miny < 0) {
        ggplt <- ggplt + geom_hline(aes_string('yintercept'=0), linetype=2, alpha=0.5)
    }
    
    r <- list(data=cor.array, plot=ggplt, separate.contexts=separate.contexts)
    return(r)
}
    
#' \code{transDist} parameter
#' 
#' Obtain an estimate for the \code{transDist} parameter (used in function \code{\link{callMethylation}}) by fitting an exponential function to the supplied correlations (from \code{\link{distanceCorrelation}}).
#' 
#' @param distcor The output produced by \code{\link{distanceCorrelation}}.
#' @param skip Skip the first n cytosines for the fitting. This can be necessary to avoid periodicity artifacts due to the context definition.
#' @param plot.parameters Whether to plot fitted parameters on to the plot or not.
#' @return A list() with fitted \code{transDist} parameters and the corresponding \code{\link[ggplot2]{ggplot}}.
#' 
#' @importFrom stats na.omit coefficients
#' @importFrom minpack.lm nlsLM
#' 
#' @export
#' @examples 
#'## Get some toy data
#'file <- system.file("data","arabidopsis_toydata.RData",
#'                     package="methimpute")
#'data <- get(load(file))
#'distcor <- distanceCorrelation(data)
#'fit <- estimateTransDist(distcor)
#'print(fit)
estimateTransDist <- function(distcor, skip=2, plot.parameters=TRUE) {

    ## Context correlation fits and plots
    contexts <- dimnames(distcor$data)[[1]]
    cor.array <- distcor$data
    maxweights <- numeric()
    params.list <- list()
    miny <- min(cor.array, na.rm = TRUE)
    dfs <- list()
    for (c1 in 1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            context.transition <- paste0(contexts[c1], '-', contexts[c2])
            if (distcor$separate.contexts) {
                if (c1 != c2) {
                    next
                }
            }
            if (c1 <= c2) {
                df <- data.frame(distance = as.numeric(dimnames(cor.array)[[3]]), correlation = cor.array[c1,c2,,'correlation'], weight = cor.array[c1,c2,,'weight'], from = contexts[c1], to = contexts[c2])
                
                ## Fit
                y <- df$correlation[(skip+1):nrow(df)]
                x <- df$distance[(skip+1):nrow(df)]
                weight <- df$weight[(skip+1):nrow(df)]
                startvalues <- list(a0 = stats::na.omit(y)[1], D = 50)
                p <- tryCatch({
                    fit <- minpack.lm::nlsLM(y ~ a0 * exp(-x/D), start=startvalues, weights=weight)
                    s <- summary(fit)
                    c <- stats::coefficients(s)
                    params <- c[1:length(startvalues)]
                    names(params) <- names(startvalues)
                    as.list(params)
                }, error = function(e) {
                    NULL
                })
                if (is.null(p)) {
                    startvalues <- list(a0 = stats::na.omit(y)[1])
                    p <- tryCatch({
                        fit <- minpack.lm::nlsLM(y ~ a0 * exp(-x/Inf), start=startvalues, weights=weight)
                        s <- summary(fit)
                        c <- stats::coefficients(s)
                        params <- c[1:length(startvalues)]
                        names(params) <- names(startvalues)
                        params <- as.list(params)
                        params$D <- Inf
                        params
                    }, error = function(e) {
                        startvalues$D <- Inf
                        startvalues
                    })
                }
                
                ## Check if we have negative D
                if (p$D <= 0) {
                    p$D <- Inf
                }
                params.list[[context.transition]] <- p
                
                ## Plot
                df$correlation.fit <- p$a0 * exp(-df$distance/p$D)
                df$logweight <- log(df$weight+1)
                dfs[[context.transition]] <- df
                maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
            }
        }
    }
    maxweight <- max(maxweights, na.rm = TRUE)
            
    ## Plot correlation
    df <- do.call(rbind, dfs)
    df$a0 <- round(sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'a0'), 2)
    df$D <- round(sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'D'), 0)
    df$params <- paste0("a0 = ", df$a0, ", D = ", df$D)
    ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
    ggplt <- ggplt + geom_line(aes_string(x='distance', y='correlation.fit'), col='blue')
    if (plot.parameters) {
        ggplt <- ggplt + geom_text(aes_string(label='params'), x=max(df$distance, na.rm = TRUE), y=max(df$correlation, na.rm = TRUE), vjust=1, hjust=1)
    }
    ggplt <- ggplt + xlab('distance in [bp]')
    ggplt <- ggplt + facet_grid(from ~ to)
    if (miny < 0) {
        ggplt <- ggplt + geom_hline(aes_string('yintercept'=0), linetype=2, alpha=0.5)
    }
    
    transDist <- sapply(params.list, '[[', 'D')
    return(list(transDist=transDist, plot=ggplt))
}
