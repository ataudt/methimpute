#' Estimate the \code{transDist} parameter
#' 
#' Estimate the \code{transDist} parameter from data.
#' 
#' @return A list() containing a data.frame and a ggplot.
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
    cors <- list()
    cor.array <- array(NA, dim=c(length(contexts), length(contexts), length(distances), 2), dimnames=list(context=contexts, context=contexts, distance=distances, what=c('correlation','weight')))
    for (i1 in distances) {
        if (i1 == distances[1]) {
            message("  for distance ", i1, appendLF=FALSE)
        } else {
            message(", ", i1, appendLF=FALSE)
        }
        ## All contexts
        ind <- which(data$distance == i1)
        ratio <- data[ind]$ratio
        ratio.shift <- data[ind-1]$ratio
        
        cor <- tryCatch(cor(ratio, ratio.shift, use='complete.obs'), error = function(e) NA, warning = function(w) NA)
        cors[[as.character(i1)]] <- data.frame(correlation = cor, weight = length(ind))
        ## Context specific
        context <- data[ind]$context
        context.shift <- data[ind-1]$context
        context.transition <- paste0(context, '-', context.shift)
        ratio.df <- data.frame(ratio, ratio.shift)
        ratio.df.splt <- split(ratio.df, context.transition)
        cor.matrix <- matrix(NA, ncol=length(contexts), nrow=length(contexts), dimnames=list(context=contexts, context=contexts))
        weight.matrix <- matrix(NA, ncol=length(contexts), nrow=length(contexts), dimnames=list(context=contexts, context=contexts))
        for (c1 in 1:length(contexts)) {
            for (c2 in c1:length(contexts)) {
                context <- paste0(contexts[c1], '-', contexts[c2])
                context.rev <- paste0(contexts[c2], '-', contexts[c1])
                df <- ratio.df.splt[[context]]
                df.rev <- ratio.df.splt[[context.rev]]
                df <- rbind(df, df.rev)
                if (!is.null(df)) {
                    cor.matrix[contexts[c2], contexts[c1]] <- tryCatch(cor(df[,1], df[,2], use='complete.obs'), error = function(e) NA, warning = function(w) NA)
                    cor.matrix[contexts[c1], contexts[c2]] <- cor.matrix[contexts[c2], contexts[c1]]
                    weight.matrix[contexts[c2], contexts[c1]] <- nrow(df)
                    weight.matrix[contexts[c1], contexts[c2]] <- weight.matrix[contexts[c2], contexts[c1]]
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
            context.transition <- paste0(contexts[c2], '-', contexts[c1])
            if (c2 > c1) {
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
    for (i1 in 1:length(ggplts)) {
        if (!is.null(ggplts[[i1]])) {
            ggplts[[i1]] <- ggplts[[i1]] + scale_alpha_continuous(name='log(weight+1)', limits=c(0,maxweight))
        }
    }
    cowplt <- cowplot::plot_grid(plotlist = ggplts, ncol=length(contexts), align='hv')
            
    ## Plot correlation
    df <- do.call(rbind, cors)
    df$distance <- distances
    df$logweight <- log(df$weight+1)
    ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
    ggplt <- ggplt + xlab('distance in [bp]')
    ggplt <- ggplt + coord_cartesian(ylim=c(0,1))
    ggplt <- ggplt + scale_alpha_continuous(name='log(weight+1)')
    
    r <- list(data=list(all = df, context = cor.array), plot=list(all = ggplt, context = cowplt))
    return(r)
}
    
#' Obtain \code{transDist} parameter.
#' 
#' Obtain an estimate for the \code{transDist} parameter in function \code{\link{callMethylation}} by fitting an exponential to the supplied correlations.
#' 
#' @param distcor The list produced by \code{\link{distanceCorrelation}}.
estimateTransDist <- function(distcor) {
  
    ## Context correlation fits and plots
    contexts <- dimnames(distcor$data$context)[[1]]
    cor.array <- distcor$data$context
    ggplts <- list()
    maxweights <- numeric()
    params.list <- list()
    for (c1 in 1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            context.transition <- paste0(contexts[c2], '-', contexts[c1])
            if (c2 > c1) {
                ggplts[context.transition] <- list(NULL)
            } else {
                df <- data.frame(distance = as.numeric(dimnames(cor.array)[[3]]), correlation = cor.array[c1,c2,,'correlation'], weight = cor.array[c1,c2,,'weight'])
                
                ## Fit
                y <- df$correlation
                x <- df$distance
                startvalues <- list(a0 = na.omit(df$correlation)[1], D = 50)
                p <- tryCatch({
                  fit <- nls(y ~ a0 * exp(-x/D), start=startvalues, weights=df$weight)
                  s <- summary(fit)
                  c <- coefficients(s)
                  params <- c[1:length(startvalues)]
                  names(params) <- names(startvalues)
                  as.list(params)
                }, error = function(e) {
                  NULL
                })
                if (is.null(p)) {
                    startvalues <- list(a0 = na.omit(df$correlation)[1])
                    p <- tryCatch({
                      fit <- nls(y ~ a0 * exp(-x/Inf), start=startvalues, weights=df$weight)
                      s <- summary(fit)
                      c <- coefficients(s)
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
                params.list[[context.transition]] <- p
                
                ## Plot
                df$correlation.fit <- p$a0 * exp(-x/p$D)
                df$logweight <- log(df$weight+1)
                ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
                ggplt <- ggplt + xlab('distance in [bp]') + ggtitle(context.transition)
                ggplt <- ggplt + coord_cartesian(ylim=c(0,1))
                ggplt <- ggplt + geom_line(aes_string(x='distance', y='correlation.fit'), col='blue')
                ggplt <- ggplt + ggtitle(paste0(context.transition, ": a0 = ", round(p$a0, 2), ", D = ", round(p$D, 2)))
                ggplts[[context.transition]] <- ggplt
                maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
            }
        }
    }
    maxweight <- max(maxweights, na.rm = TRUE)
    for (i1 in 1:length(ggplts)) {
        if (!is.null(ggplts[[i1]])) {
            ggplts[[i1]] <- ggplts[[i1]] + scale_alpha_continuous(name='log(weight+1)', limits=c(0,maxweight))
        }
    }
    cowplt <- cowplot::plot_grid(plotlist = ggplts, ncol=length(contexts), align='hv')
            
    return(list(params=params.list, plot=cowplt))
}