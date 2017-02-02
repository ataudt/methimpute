#' Distance correlation
#' 
#' Compute the distance correlation from a \code{\link{methimputeData}} object.
#' 
#' @param data A \code{\link{methimputeData}} object.
#' @param distances An integer vector specifying the distances for which the correlation will be calculated.
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
distanceCorrelation <- function(data, distances=0:50) {
    
    ## Contexts
    contexts <- intersect(levels(data$context), unique(data$context))
    
    ## Add ratio column
    data$ratio <- data$counts[,'methylated'] / data$counts[,'total']
    ### Add distance to bins ###
    ptm <- startTimedMessage("Adding distance ...")
    data$distance <- c(-1, start(data)[-1] - end(data)[-length(data)] - 1)
    data$distance[data$distance < 0] <- Inf 
    stopTimedMessage(ptm)
  
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
        ratio <- data[ind]$ratio
        ratio.shift <- data[ind-1]$ratio
        # ## All contexts
        # cor <- tryCatch(cor(ratio, ratio.shift, use='complete.obs'), error = function(e) NA, warning = function(w) NA)
        # cors[[as.character(i1)]] <- data.frame(correlation = cor, weight = length(ind))
        ## Context specific
        context <- data$context[ind]
        context.shift <- data$context[ind-1]
        context.transition <- paste0(context, '-', context.shift)
        cor.matrix <- matrix(NA, ncol=length(contexts), nrow=length(contexts), dimnames=list(context=contexts, context=contexts))
        weight.matrix <- matrix(NA, ncol=length(contexts), nrow=length(contexts), dimnames=list(context=contexts, context=contexts))
        for (c1 in 1:length(contexts)) {
            for (c2 in c1:length(contexts)) {
                context <- paste0(contexts[c1], '-', contexts[c2])
                context.rev <- paste0(contexts[c2], '-', contexts[c1])
                mask <- context.transition == context | context.transition == context.rev
                ratio1 <- ratio[mask]
                ratio2 <- ratio.shift[mask]
                if (!is.null(df)) {
                    cor.matrix[contexts[c1], contexts[c2]] <- tryCatch(cor(ratio1, ratio2, use='complete.obs'), error = function(e) NA, warning = function(w) NA)
                    cor.matrix[contexts[c2], contexts[c1]] <- cor.matrix[contexts[c1], contexts[c2]]
                    weight.matrix[contexts[c1], contexts[c2]] <- length(ratio1)
                    weight.matrix[contexts[c2], contexts[c1]] <- weight.matrix[contexts[c1], contexts[c2]]
                }
            }
        }
        cor.array[,,as.character(i1),'correlation'] <- cor.matrix
        cor.array[,,as.character(i1),'weight'] <- weight.matrix
    }
    message("\nFinished calculating correlations in", appendLF = FALSE)
    stopTimedMessage(ptm)
    
    ## Context correlation plots
    ggplts <- list()
    maxweights <- numeric()
    for (c1 in 1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            context.transition <- paste0(contexts[c1], '-', contexts[c2])
            if (c1 > c2) {
                ggplts[context.transition] <- list(NULL)
            } else {
                df <- data.frame(distance = distances, correlation = cor.array[c1,c2,,'correlation'], weight = cor.array[c1,c2,,'weight'])
                df$logweight <- log(df$weight+1)
                ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
                ggplt <- ggplt + xlab('distance in [bp]') + ggtitle(context.transition)
                ggplt <- ggplt + coord_cartesian(ylim=c(0,1))
                ggplts[[context.transition]] <- ggplt
                maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
            }
        }
    }
    maxweight <- max(maxweights, na.rm = TRUE)
    miny <- min(cor.array, na.rm = TRUE)
    for (i1 in 1:length(ggplts)) {
        if (!is.null(ggplts[[i1]])) {
            ggplts[[i1]] <- ggplts[[i1]] + scale_alpha_continuous(name='log(weight+1)', limits=c(0,maxweight))
            ggplts[[i1]] <- ggplts[[i1]] + coord_cartesian(ylim=c(min(0, miny), 1))
        }
    }
    cowplt <- suppressWarnings( cowplot::plot_grid(plotlist = ggplts, ncol=length(contexts), align='hv') )
            
    # ## Plot correlation
    # df <- do.call(rbind, cors)
    # df$distance <- distances
    # df$logweight <- log(df$weight+1)
    # ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
    # ggplt <- ggplt + xlab('distance in [bp]')
    # ggplt <- ggplt + coord_cartesian(ylim=c(0,1))
    # ggplt <- ggplt + scale_alpha_continuous(name='log(weight+1)')
    # r <- list(data=list(all = df, context = cor.array), plot=list(all = ggplt, context = cowplt))
    
    r <- list(data=cor.array, plot=cowplt)
    return(r)
}
    
#' \code{transDist} parameter
#' 
#' Obtain an estimate for the \code{transDist} parameter (used in function \code{\link{callMethylation}}) by fitting an exponential function to the supplied correlations (from \code{\link{distanceCorrelation}}).
#' 
#' @param distcor The output produced by \code{\link{distanceCorrelation}}.
#' @param skip Skip the first n cytosines for the fitting. This can be necessary to avoid periodicity artifacts due to the context definition.
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
estimateTransDist <- function(distcor, skip=2) {
  
    ## Context correlation fits and plots
    contexts <- dimnames(distcor$data)[[1]]
    cor.array <- distcor$data
    ggplts <- list()
    maxweights <- numeric()
    params.list <- list()
    miny <- min(cor.array, na.rm = TRUE)
    for (c1 in 1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            context.transition <- paste0(contexts[c1], '-', contexts[c2])
            if (c1 > c2) {
                ggplts[context.transition] <- list(NULL)
            } else {
                df <- data.frame(distance = as.numeric(dimnames(cor.array)[[3]]), correlation = cor.array[c1,c2,,'correlation'], weight = cor.array[c1,c2,,'weight'])
                
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
                ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
                ggplt <- ggplt + xlab('distance in [bp]') + ggtitle(context.transition)
                ggplt <- ggplt + geom_line(aes_string(x='distance', y='correlation.fit'), col='blue')
                ggplt <- ggplt + ggtitle(paste0(context.transition, ": a0 = ", round(p$a0, 2), ", D = ", round(p$D)))
                ggplts[[context.transition]] <- ggplt
                maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
            }
        }
    }
    maxweight <- max(maxweights, na.rm = TRUE)
    for (i1 in 1:length(ggplts)) {
        if (!is.null(ggplts[[i1]])) {
            ggplts[[i1]] <- ggplts[[i1]] + scale_alpha_continuous(name='log(weight+1)', limits=c(0,maxweight))
            ggplts[[i1]] <- ggplts[[i1]] + coord_cartesian(ylim=c(min(0, miny), 1))
        }
    }
    cowplt <- suppressWarnings( cowplot::plot_grid(plotlist = ggplts, ncol=length(contexts), align='hv') )
            
    transDist <- sapply(params.list, '[[', 'D')
    return(list(transDist=transDist, plot=cowplt))
}
