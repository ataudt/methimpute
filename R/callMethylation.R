#' Call methylation status
#' 
#' Call methylation status of cytosines (or bins) with a Hidden Markov Model.
#' 
#' The Hidden Markov model uses a binomial test for the emission densities. Transition probabilities are modeled with a distance dependent decay, specified by the parameter \code{transDist}.
#' 
#' @param data A \code{\link{methimputeData}} object.
#' @param fit.on.chrom A character vector specifying the chromosomes on which the HMM will be fitted.
#' @param transDist The decaying constant for the distance-dependent transition matrix. Either a single numeric or a named numeric vector, where the vector names correspond to the transition contexts. Such a vector can be obtained from \code{\link{estimateTransDist}}.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param max.time Maximum running time in seconds for the Baum-Welch algorithm.
#' @param max.iter Maximum number of iterations for the Baum-Welch algorithm.
#' @param count.cutoff A cutoff for the counts to remove artificially high counts from mapping artifacts. Set to \code{Inf} to disable this filtering (not recommended).
#' @param verbosity An integer from 1 to 5 specifying the verbosity of the fitting procedure. Values > 1 are only for debugging.
#' @param num.threads Number of CPU to use for the computation. Parallelization is implemented on the number of states, which is 2 or 3, so setting \code{num.threads > 3} will not give additional performance increase.
#' @param initial.params A \code{\link{methimputeBinomialHMM}} object. This parameter is useful to continue the fitting procedure for a \code{\link{methimputeBinomialHMM}} object.
#' @param include.intermediate A logical specifying wheter or not the intermediate component should be included in the HMM.
#' @param update One of \code{c("independent", "constrained")}. If \code{update="independent"} probability parameters for the binomial test will be updated independently. If \code{update="constrained"} the probability parameter of the intermediate component will be constrained to the mean of the unmethylated and the methylated component.
#' @param min.reads The minimum number of reads that a position must have to contribute in the Baum-Welch fitting procedure.
#' @return A \code{\link{methimputeBinomialHMM}} object.
#' 
#' @export
#' @examples
#'## Get some toy data
#'file <- system.file("data","arabidopsis_toydata.RData", package="methimpute")
#'data <- get(load(file))
#'print(data)
#'model <- callMethylation(data)
#'print(model)
callMethylation <- function(data, fit.on.chrom=NULL, transDist=Inf, eps=1, max.time=Inf, max.iter=Inf, count.cutoff=500, verbosity=1, num.threads=2+include.intermediate, initial.params=NULL, include.intermediate=FALSE, update='independent', min.reads=0) {
  
    ### Input checks ###
    if (!is.null(fit.on.chrom)) {
        if (!fit.on.chrom %in% seqlevels(data)) {
            stop("Cannot find 'fit.on.chrom' = ", fit.on.chrom, " in the data.")
        }
    }
  
    ### Add distance and transition context to bins ###
    data$distance <- addDistance(data)
    data$transitionContext <- addTransitionContext(data)
    
    ### Assign variables ###
    update <- factor(update, levels=c('independent', 'constrained'))
    if (is.na(update)) { stop("Argument 'update' must be one of c('independent', 'constrained').") }
    contexts <- intersect(levels(data$context), unique(data$context))
    ncontexts <- length(contexts)
    transitionContexts <- levels(data$transitionContext)
    transDistvec <- rep(Inf, length(transitionContexts))
    names(transDistvec) <- transitionContexts
    if (is.null(names(transDist))) {
        transDistvec[1:length(transDistvec)] <- transDist
    } else {
        if (any(! names(transDist) %in% names(transDistvec))) {
            stop("Names of 'transDist' must be ", paste0(names(transDistvec), collapse=', '), ".")
        }
        transDistvec[names(transDist)] <- transDist
        rev.names <- sapply(strsplit(names(transDist), '-'), function(x) { paste0(rev(x), collapse = '-')})
        transDistvec[rev.names] <- transDist
    }
    if (!include.intermediate) {
        states <- c('Unmethylated', 'Methylated')
    } else {
        states <- c('Unmethylated', 'Intermediate', 'Methylated')
    }
    numstates <- length(states)
    counts <- data$counts
    distances <- data$distance
    context <- factor(data$context, levels=contexts)
    transitionContext <- factor(data$transitionContext, levels=transitionContexts)
    
    
    ## Filter counts by cutoff
    mask <- counts[,'total'] > count.cutoff
    counts[mask,] <- round(sweep(x = counts[mask,, drop=FALSE], MARGIN = 1, STATS = counts[mask,'total', drop=FALSE]/count.cutoff, FUN = '/'))
    data$counts <- counts # assign it now to have filtered values in there
    
    ## Subset by chromosomes
    if (!is.null(fit.on.chrom)) {
        mask <- as.logical(data@seqnames %in% fit.on.chrom)
        counts <- counts[mask,]
        context <- context[mask]
        transitionContext <- transitionContext[mask]
        distances <- distances[mask]
    }
  
    ### Initial probabilities ###
    if (is.null(initial.params)) {
        transProbs_initial <- list()
        for (context.transition in transitionContexts) {
            s <- 0.9
            transProbs_initial[[context.transition]] <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=states, to=states))
            diag(transProbs_initial[[context.transition]]) <- s
        }
        startProbs_initial <- rep(1/numstates, numstates)
        names(startProbs_initial) <- states
        
        ### Initialization of emission distributions ###
        ep <- list()
        probUN.start <- 0.01
        probM.start <- 0.9
        if (!include.intermediate) {
            probs <- rep(probUN.start, ncontexts)
            names(probs) <- contexts
            ep[[states[1]]] <- data.frame(prob=probs)
            probs <- rep(probM.start, ncontexts)
            names(probs) <- contexts
            ep[[states[2]]] <- data.frame(prob=probs)
        } else {
            probs <- rep(probUN.start, ncontexts)
            names(probs) <- contexts
            ep[[states[1]]] <- data.frame(prob=probs)
            probs <- rep(0.5*(probUN.start+probM.start), ncontexts)
            names(probs) <- contexts
            ep[[states[2]]] <- data.frame(prob=probs)
            probs <- rep(probM.start, ncontexts)
            names(probs) <- contexts
            ep[[states[3]]] <- data.frame(prob=probs)
        }
        names(ep) <- states
        emissionParams_initial <- ep
    } else {
        model.init <- loadFromFiles(initial.params, check.class='methimputeBinomialHMM')[[1]]
        transProbs_initial <- model.init$params$transProbs
        startProbs_initial <- model.init$params$startProbs
        emissionParams_initial <- model.init$params$emissionParams
    }
  
    ### Define parameters for C function ###
    params <- list()
    params$startProbs_initial <- startProbs_initial
    params$transProbs_initial <- transProbs_initial
    params$emissionParams_initial <- emissionParams_initial
    params$transDist <- transDistvec
    params$eps <- eps
    params$maxtime <- max.time
    params$maxiter <- max.iter
    params$minreads <- min.reads
    params$verbosity <- verbosity
    params$numThreads <- num.threads
    
    ### Call the C function ###
    on.exit(cleanup())
    if (is.null(fit.on.chrom)) {
        ptm <- startTimedMessage("Baum-Welch: Fitting HMM parameters\n")
        hmm <- fitBinomialTestHMMcontextTransition(counts_total=counts[,'total'], counts_meth=counts[,'methylated'], context=as.integer(context)-1, transitionContext=as.integer(transitionContext)-1, distances=distances, params=params, algorithm=1, update_procedure=update)
        message("Time spent in Baum-Welch:", appendLF=FALSE)
        stopTimedMessage(ptm)
        if (hmm$error == "") {
            ## Cast convergence info
            parray <- array(NA, dim=c(length(contexts), length(hmm$convergenceInfo$logliks), length(states)), dimnames=list(context=contexts, iteration=0:(length(hmm$convergenceInfo$logliks)-1), status=states))
            parray[,,'Unmethylated'] <- hmm$convergenceInfo$parameterInfo$probsUN
            parray[,,'Methylated'] <- hmm$convergenceInfo$parameterInfo$probsM
            if ('Intermediate' %in% states) {
                parray[,,'Intermediate'] <- (parray[,,'Unmethylated'] + parray[,,'Methylated']) / 2
            }
            convergenceInfo <- hmm$convergenceInfo
            convergenceInfo$parameterInfo <- parray[,-1,]
            convergenceInfo$logliks <- convergenceInfo$logliks[-1]
            convergenceInfo$dlogliks <- convergenceInfo$dlogliks[-1]
        }
    } else {
        ptm <- startTimedMessage("Baum-Welch: Fitting HMM parameters\n")
        message(" ... on chromosomes ", paste0(fit.on.chrom, collapse=', '))
        hmm <- fitBinomialTestHMMcontextTransition(counts_total=counts[,'total'], counts_meth=counts[,'methylated'], context=as.integer(context)-1, transitionContext=as.integer(transitionContext)-1, distances=distances, params=params, algorithm=1, update_procedure=update)
        message("Time spent in Baum-Welch:", appendLF=FALSE)
        stopTimedMessage(ptm)
        if (hmm$error == "") {
            ## Cast convergence info
            convergenceInfo <- list(prefit=list(), fit=list())
            convergenceInfo$prefit$fit.on.chrom <- fit.on.chrom
            parray <- array(NA, dim=c(length(contexts), length(hmm$convergenceInfo$logliks), length(states)), dimnames=list(context=contexts, iteration=0:(length(hmm$convergenceInfo$logliks)-1), status=states))
            parray[,,'Unmethylated'] <- hmm$convergenceInfo$parameterInfo$probsUN
            parray[,,'Methylated'] <- hmm$convergenceInfo$parameterInfo$probsM
            if ('Intermediate' %in% states) {
                parray[,,'Intermediate'] <- (parray[,,'Unmethylated'] + parray[,,'Methylated']) / 2
            }
            convergenceInfo$prefit <- hmm$convergenceInfo
            convergenceInfo$prefit$parameterInfo <- parray[,-1,]
            convergenceInfo$prefit$logliks <- convergenceInfo$prefit$logliks[-1]
            convergenceInfo$prefit$dlogliks <- convergenceInfo$prefit$dlogliks[-1]
            ## Redo for all chromosomes
            counts <- data$counts
            context <- data$context
            transitionContext <- data$transitionContext
            distances <- data$distance
            params2 <- list()
            params2$startProbs_initial <- hmm$startProbs
            params2$transProbs_initial <- hmm$transProbs
            params2$emissionParams_initial <- hmm$emissionParams
            params2$transDist <- transDistvec
            params2$eps <- eps
            params2$maxtime <- max.time
            params2$maxiter <- max.iter
            params2$minreads <- min.reads
            params2$verbosity <- verbosity
            params2$numThreads <- num.threads
            ptm <- startTimedMessage("Forward-Backward: Obtaining state sequence - no updates\n")
            message(" ... on all chromosomes")
            hmm <- fitBinomialTestHMMcontextTransition(counts_total=counts[,'total'], counts_meth=counts[,'methylated'], context=as.integer(context)-1, transitionContext=as.integer(transitionContext)-1, distances=distances, params=params2, algorithm=2, update_procedure=update)
            message("Time spent in Forward-Backward:", appendLF=FALSE)
            stopTimedMessage(ptm)
            ## Cast convergence info
            convergenceInfo$fit <- hmm$convergenceInfo
        }
    }
    
    ### Construct result object ###
    ptm <- startTimedMessage("Compiling results ...")
    r <- list()
    class(r) <- "methimputeBinomialHMM"
    if (hmm$error == "") {
        r$convergenceInfo <- convergenceInfo
        names(hmm$weights) <- states
        r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, transDist=hmm$transDist, emissionParams=hmm$emissionParams, count.cutoff=count.cutoff)
        r$params.initial <- params
        # States and posteriors
        rownames(hmm$posteriors) <- states
        data$posteriorMax <- NA
        for (i1 in 0:(nrow(hmm$posteriors)-1)) {
            mask <- hmm$states == i1
            data$posteriorMax[mask] <- hmm$posteriors[i1+1,mask]
        }
        data$posteriorMeth <- hmm$posteriors["Methylated",]
        data$posteriorUnmeth <- hmm$posteriors["Unmethylated",]
        data$status <- factor(states, levels=states)[hmm$states+1]
        # Recalibrated methylation levels
        if (include.intermediate) {
            data$rc.meth.lvl <- r$params$emissionParams$Unmethylated[data$context,] * data$posteriorUnmeth + r$params$emissionParams$Intermediate[data$context,] * (1 - data$posteriorMeth - data$posteriorUnmeth) + r$params$emissionParams$Methylated[data$context,] * data$posteriorMeth
        } else {
            data$rc.meth.lvl <- r$params$emissionParams$Unmethylated[data$context,] * data$posteriorUnmeth + r$params$emissionParams$Methylated[data$context,] * data$posteriorMeth
        }
        ## Recalibrated counts
        data$rc.counts <- data$counts
        # Transform posteriorMax to uniform
        ecdf.postmax <- ecdf(data$posteriorMax)
        upostmax <- ecdf.postmax(data$posteriorMax)
        # Transform to count distribution
        x <- seq(0, 1, by=1e-5)
        q <- quantile(data$counts[,2], probs = x)
        inverse.ecdf <- stepfun(x = x, y = c(q, max(q)))
        cpostmax <- inverse.ecdf(upostmax)
        # Assign to counts
        data$rc.counts[,2] <- round(cpostmax)
        data$rc.counts[,1] <- round(data$rc.counts[,2] * data$rc.meth.lvl)
        ## Segmentation
        data.collapse <- data
        mcols(data.collapse) <- NULL
        data.collapse$status <- data$status
        df <- suppressMessages( collapseBins(as.data.frame(data.collapse), column2collapseBy = 'status') )
        segments <- methods::as(df, 'GRanges')
        seqlengths(segments) <- seqlengths(data)[seqlevels(segments)]
        ## Context-dependent weights
        r$params$weights <- list()
        for (context in contexts) {
            r$params$weights[[context]] <- rowMeans(hmm$posteriors[,data$context==context])
        }
    } else {
        warning("Baum-Welch aborted: ", hmm$error)
    }
    r$data <- data
    r$segments <- segments
    stopTimedMessage(ptm)
    
    return(r)
}


#' Call methylation status
#' 
#' Call methylation status of cytosines (or bins) with a binomial test.
#' 
#' The function uses a binomial test with the specified \code{conversion.rate}. P-values are then multiple testing corrected with the Benjamini & Yekutieli procedure. Methylated positions are selected with the \code{p.threshold}.
#' 
#' @param data A \code{\link{methimputeData}} object.
#' @param conversion.rate A conversion rate between 0 and 1.
#' @param min.coverage Minimum coverage to consider for the binomial test.
#' @param p.threshold Significance threshold between 0 and 1.
#' @return A vector with methylation statuses.
#' 
#' @examples
#'## Get some toy data
#'file <- system.file("data","arabidopsis_toydata.RData", package="methimpute")
#'data <- get(load(file))
#'data$binomial <- binomialTestMethylation(data, conversion.rate=0.998)
#'
binomialTestMethylation <- function(data, conversion.rate, min.coverage=3, p.threshold=0.05) {
  
    p <- pbinom(q = data$counts[,'methylated']-1, size = data$counts[,'total'], prob = conversion.rate, lower.tail = FALSE)
    p[data$counts[,'total'] < min.coverage] <- NA
    p <- p.adjust(p, method = 'BY')
    levels <- c("Unmethylated", "Methylated")
    methylated <- factor(levels[c(1,2)][(p <= p.threshold)+1], levels=levels)
    return(methylated)
  
}

#' #' Make a multivariate segmentation
#' #' 
#' #' Make a multivariate segmentation by fitting the transition probabilities of a multivariate Hidden Markov Model.
#' #' 
#' #' The function will use the provided univariate Hidden Markov Models (HMMs) to build the multivariate emission density. Number of states are also taken from the combination of univaraite HMMs.
#' #' 
#' #' @return A list with fitted parameters, posteriors.
#' multivariateSegmentation <- function(models, ID, fit.on.chrom=NULL, transDist=700, eps=0.01, max.time=Inf, max.iter=Inf, max.states=Inf, verbosity=1, num.threads=1) {
#'   
#'     ### Input checks ###
#'     if (is.null(max.time)) {
#'         max.time <- Inf
#'     }
#'     if (is.null(max.iter)) {
#'         max.iter <- Inf
#'     }
#'     ID.check <- ID
#'   
#'     ### Assign variables ###
#'     models <- loadFromFiles(models, check.class = 'NcomponentHMM')
#'     components <- lapply(models, function(model) { levels(model$data$state) })
#'     statedef <- permute(components)
#'     # Data object
#'     data <- models[[1]]$data
#'     mcols(data) <- NULL
#'     data$distance <- models[[1]]$data$distance
#'     
#'     ### Prepare the multivariate ###
#'     messageU("Multivariate HMM: preparing fitting procedure", underline="=", overline="=")
#'     cormat <- prepareMultivariate(data=models[[1]]$data, hmms=models, statedef=statedef, max.states=max.states)
#'     data$observable <- cormat$observable
#'     statenames <- cormat$mapping
#'     statenames <- lapply(strsplit(statenames, ' '), function(x) { paste(paste0(names(models), "=", x), collapse=" ") })
#'     numstates <- length(statenames)
#'     statedef <- cormat$statedef
#'     
#'     ### Initial probabilities ###
#'     s <- 0.9
#'     transProbs_initial <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=statenames, to=statenames))
#'     diag(transProbs_initial) <- s
#'     startProbs_initial <- rep(1/numstates, numstates)
#'     names(startProbs_initial) <- statenames
#'     
#'     ### Define parameters for C function ###
#'     params <- list()
#'     params$startProbs_initial <- startProbs_initial
#'     params$transProbs_initial <- transProbs_initial
#'     params$emissionParamsList <- cormat$emissionParams
#'     params$transDist <- transDist
#'     params$eps <- eps
#'     params$maxtime <- max.time
#'     params$maxiter <- max.iter
#'     params$verbosity <- verbosity
#'     params$numThreads <- num.threads
#'     params$correlationMatrixInverse <- cormat$correlationMatrixInverse
#'     params$determinant <- cormat$determinant
#'     params$statedef <- cormat$statedef
#'     
#'     ### Fit the HMM ###
#'     message("Baum-Welch: Fitting parameters for multivariate HMM")
#'     on.exit(cleanup())
#'     hmm <- fitMultiHMM(data$observable, data$distance, params)
#'     
#'     ### Construct result object ###
#'     ptm <- startTimedMessage("Compiling results ...")
#'     r <- list()
#'     class(r) <- 'NcomponentMultiHMM'
#'     r$ID <- ID
#'     if (hmm$error == "") {
#'         r$convergenceInfo <- hmm$convergenceInfo
#'         names(hmm$weights) <- statenames
#'         r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, transDist=hmm$transDist, emissionParamsList=params$emissionParamsList, weights=hmm$weights)
#'         r$params.initial <- params
#'         # States and posteriors
#'         data$posteriorMax <- NA
#'         for (i1 in 0:(nrow(hmm$posteriors)-1)) {
#'             mask <- hmm$states == i1
#'             data$posteriorMax[mask] <- hmm$posteriors[i1+1,mask]
#'         }
#'         rownames(hmm$posteriors) <- statenames
#'         data$posteriors <- t(hmm$posteriors)
#'         data$state <- factor(statenames, levels=statenames)[hmm$states+1]
#'         ## Make segmentation
#'         df <- as.data.frame(data)
#'         df <- df[, c(names(df)[1:5], 'state')]
#'         segments <- suppressMessages( collapseBins(df, column2collapseBy = 'state') )
#'         segments <- methods::as(segments, 'GRanges')
#'         seqlevels(segments) <- seqlevels(data)
#'         seqlengths(segments) <- seqlengths(data)[seqlevels(segments)]
#'     }
#'     r$segments <- segments
#'     r$data <- data
#'     stopTimedMessage(ptm)
#'     
#'     return(r)
#' }
#' 

#' #' Fit a two-component HMM
#' #' 
#' #' Fit a two-component Hidden Markov Model to the supplied counts. The transition matrix is distance-dependent with exponential decaying constant \code{transDist}. Components are modeled as negative binomial distributions.
#' #' 
#' #' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata column 'counts' (or any other column specified as \code{observable}).
#' #' @param observable A character naming the metadata column of \code{data} that will serve as observable for the HMM.
#' #' @param fit.on.chrom A character vector giving the chromosomes on which the HMM will be fitted.
#' #' @param transDist The exponential decaying constant for the distance-dependent transition matrix. Should be given in the same units as \code{distances}.
#' #' @param eps Convergence threshold for the Baum-Welch algorithm.
#' #' @param max.time Maximum running time in seconds for the Baum-Welch algorithm.
#' #' @param max.iter Maximum number of iterations for the Baum-Welch algorithm.
#' #' @param count.cutoff A cutoff for the \code{observable}. Lower cutoffs will speed up the fitting procedure and improve convergence in some cases. Set to \code{Inf} to disable this filtering.
#' #' @param verbosity Integer from c(0,1) specifying the verbosity of the fitting procedure.
#' #' @param num.threads Number of CPU to use for the computation.
#' #' @return A list with fitted parameters, posteriors, and the input parameters.
#' #' 
#' fitSignalBackground <- function(data, observable='counts', fit.on.chrom=NULL, transDist=700, eps=0.01, max.time=Inf, max.iter=Inf, cutoff=1000, verbosity=1, num.threads=1) {
#'   
#'     ### Input checks ###
#'     if (is.null(max.time)) {
#'         max.time <- Inf
#'     }
#'     if (is.null(max.iter)) {
#'         max.iter <- Inf
#'     }
#'   
#'     ### Add distance to bins ###
#'     data$distance <- addDistance(data)
#'     
#'     ### Assign variables ###
#'     states <- c("background", "signal")
#'     numstates <- length(states)
#'     counts <- mcols(data)[,observable]
#'     counts <- as.matrix(counts)
#'     distances <- data$distance
#'     
#'     ## Filter counts by cutoff
#'     mask <- counts[,'total'] > count.cutoff
#'     counts[mask,] <- round(sweep(x = counts[mask,], MARGIN = 1, STATS = counts[mask,'total']/count.cutoff, FUN = '/'))
#'     data$observable <- counts # assign it now to have filtered values in there
#'     colnames(data$observable) <- observable
#'     
#'     ## Subset by chromosomes
#'     if (!is.null(fit.on.chrom)) {
#'         counts <- counts[as.logical(data@seqnames %in% fit.on.chrom)]
#'     }
#'   
#'     ### Initial probabilities ###
#'     s <- 0.9
#'     transProbs_initial <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=states, to=states))
#'     diag(transProbs_initial) <- s
#'     startProbs_initial <- rep(1/numstates, numstates)
#'     names(startProbs_initial) <- states
#'     
#'     ### Initialization of emission distributions ###
#'     counts.greater.0 <- counts[counts>0]
#'     mean.counts <- mean(counts.greater.0)
#'     var.counts <- var(counts.greater.0)
#'     ep <- list()
#'     ep[["background"]] <- data.frame(type='dnbinom', mu=mean.counts, var=var.counts)
#'     ep[["signal"]] <- data.frame(type='dnbinom', mu=mean.counts*1.5, var=var.counts*2)
#'     ep <- do.call(rbind, ep)
#'     # Make sure variance is bigger than mean for negative binomial
#'     ep$var[ep$mu >= ep$var] <- ep$mu[ep$mu >= ep$var] + 1
#'     ep$size <- dnbinom.size(ep$mu, ep$var)
#'     ep$prob <- dnbinom.prob(ep$mu, ep$var)
#'     emissionParams_initial <- ep
#'   
#'     ### Define parameters for C function ###
#'     params <- list()
#'     params$startProbs_initial <- startProbs_initial
#'     params$transProbs_initial <- transProbs_initial
#'     params$emissionParams_initial <- emissionParams_initial
#'     params$transDist <- transDist
#'     params$eps <- eps
#'     params$maxtime <- max.time
#'     params$maxiter <- max.iter
#'     params$verbosity <- verbosity
#'     params$numThreads <- num.threads
#'     
#'     ### Call the C function ###
#'     on.exit(cleanup())
#'     if (is.null(fit.on.chrom)) {
#'         ptm <- startTimedMessage("Baum-Welch: Fitting HMM parameters\n")
#'         hmm <- fitHMM(counts=counts, distances=distances, params=params, algorithm=1)
#'         message("Time spent in Baum-Welch:", appendLF=FALSE)
#'         stopTimedMessage(ptm)
#'     } else {
#'         ptm <- startTimedMessage("Baum-Welch: Fitting HMM parameters\n")
#'         message(" ... on chromosomes ", paste0(fit.on.chrom, collapse=', '))
#'         hmm <- fitHMM(counts=counts, distances=distances, params=params, algorithm=1)
#'         message("Time spent in Baum-Welch:", appendLF=FALSE)
#'         stopTimedMessage(ptm)
#'         counts <- as.integer(data$observable)
#'         params2 <- list()
#'         params2$startProbs_initial <- hmm$startProbs
#'         params2$transProbs_initial <- hmm$transProbs
#'         params2$emissionParams_initial <- hmm$emissionParams
#'         params2$transDist <- transDist
#'         params2$eps <- eps
#'         params2$maxtime <- max.time
#'         params2$maxiter <- max.iter
#'         params2$verbosity <- verbosity
#'         params2$numThreads <- num.threads
#'         ptm <- startTimedMessage("Viterbi: Obtaining state sequence\n")
#'         message(" ... on all chromosomes")
#'         hmm <- fitHMM(counts=counts, distances=distances, params=params2, algorithm=2)
#'         message("Time spent in Viterbi:", appendLF=FALSE)
#'         stopTimedMessage(ptm)
#'     }
#'     
#'     ### Construct result object ###
#'     ptm <- startTimedMessage("Compiling results ...")
#'     r <- list()
#'     if (hmm$error == "") {
#'         r$convergenceInfo <- hmm$convergenceInfo
#'         names(hmm$weights) <- states
#'         r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, transDist=hmm$transDist, emissionParams=hmm$emissionParams, weights=hmm$weights)
#'         r$params.initial <- params
#'         # States and posteriors
#'         data$posteriorMax <- NA
#'         for (i1 in 0:(nrow(hmm$posteriors)-1)) {
#'             mask <- hmm$states == i1
#'             data$posteriorMax[mask] <- hmm$posteriors[i1+1,mask]
#'         }
#'         rownames(hmm$posteriors) <- states
#'         data$posteriors <- t(hmm$posteriors)
#'         data$state <- factor(states, levels=states)[hmm$states+1]
#'     }
#'     r$data <- data
#'     stopTimedMessage(ptm)
#'     
#'     return(r)
#' }
#' 
#' 
#' 
#' 
#' #' Fit an n-component HMM
#' #' 
#' #' Fit an n-component Hidden Markov Model to the supplied counts. The transition matrix is distance-dependent with exponential decaying constant \code{transDist} (only relevant in non-consecutive bins). The zero-th component is a delta distribution to account for empty bins, and all other n-components are modeled as negative binomial distributions.
#' #' 
#' #' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'distance' and 'counts' (or any other column specified as \code{observable}).
#' #' @param states An integer vector giving the states for the Hidden Markov Model. State '0' will be modeled by a delta distribution, all other states ('1','2','3',...) with negative binomial distributions.
#' #' @param observable A character naming the column of \code{data} that will serve as observable for the HMM.
#' #' @param fit.on.chrom A character vector giving the chromosomes on which the HMM will be fitted.
#' #' @param transDist The exponential decaying constant for the distance-dependent transition matrix. Should be given in the same units as \code{distances}.
#' #' @param eps Convergence threshold for the Baum-Welch algorithm.
#' #' @param max.time Maximum running time in seconds for the Baum-Welch algorithm.
#' #' @param max.iter Maximum number of iterations for the Baum-Welch algorithm.
#' #' @param count.cutoff A cutoff for the \code{observable}. Lower cutoffs will speed up the fitting procedure and improve convergence in some cases. Set to \code{Inf} to disable this filtering.
#' #' @param verbosity Integer from c(0,1) specifying the verbosity of the fitting procedure.
#' #' @param num.threads Number of CPU to use for the computation.
#' #' @param initial.params An \code{\link{NcomponentHMM}} or a file that contains such an object. Parameters from this model will be used for initialization of the fitting procedure.
#' #' @return A list with fitted parameters, posteriors, and the input parameters.
#' #' 
#' fitNComponentHMM <- function(data, states=0:5, observable='counts', fit.on.chrom=NULL, transDist=700, eps=0.01, max.time=Inf, max.iter=Inf, count.cutoff=500, verbosity=1, num.threads=1, initial.params=NULL) {
#'   
#'     ### Input checks ###
#'     if (is.null(max.time)) {
#'         max.time <- Inf
#'     }
#'     if (is.null(max.iter)) {
#'         max.iter <- Inf
#'     }
#'   
#'     ### Add distance to bins ###
#'     data$distance <- addDistance(data)
#'     
#'     ### Assign variables ###
#'     numstates <- length(states)
#'     counts <- mcols(data)[,observable]
#'     counts <- as.matrix(counts)
#'     distances <- data$distance
#'     
#'     ## Filter counts by cutoff
#'     mask <- counts[,'total'] > count.cutoff
#'     counts[mask,] <- round(sweep(x = counts[mask,], MARGIN = 1, STATS = counts[mask,'total']/count.cutoff, FUN = '/'))
#'     data$observable <- counts # assign it now to have filtered values in there
#'     colnames(data$observable) <- observable
#'     
#'     ## Subset by chromosomes
#'     if (!is.null(fit.on.chrom)) {
#'         counts <- counts[as.logical(data@seqnames %in% fit.on.chrom)]
#'     }
#'   
#'     ### Initial probabilities ###
#'     if (is.null(initial.params)) {
#'         s <- 0.9
#'         transProbs_initial <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=states, to=states))
#'         diag(transProbs_initial) <- s
#'         startProbs_initial <- rep(1/numstates, numstates)
#'         names(startProbs_initial) <- states
#'         
#'         ### Initialization of emission distributions ###
#'         counts.greater.0 <- counts[counts>0]
#'         mean.counts <- mean(counts.greater.0)
#'         var.counts <- var(counts.greater.0)
#'         ep <- list()
#'         for (state in setdiff(states,0)) {
#'             ep[[as.character(state)]] <- data.frame(type='dnbinom', mu=mean.counts*state/2, var=var.counts*state/2)
#'         }
#'         ep <- do.call(rbind, ep)
#'         # Make sure variance is bigger than mean for negative binomial
#'         ep$var[ep$mu >= ep$var] <- ep$mu[ep$mu >= ep$var] + 1
#'         ep$size <- dnbinom.size(ep$mu, ep$var)
#'         ep$prob <- dnbinom.prob(ep$mu, ep$var)
#'         ## Add zero-inflation
#'         if (0 %in% states) {
#'             ep <- rbind(data.frame(type='delta', mu=0, var=0, size=NA, prob=NA), ep)
#'         }
#'         rownames(ep) <- states
#'         emissionParams_initial <- ep
#'     } else {
#'         model.init <- loadFromFiles(initial.params, check.class='NcomponentHMM')[[1]]
#'         transProbs_initial <- model.init$params$transProbs
#'         startProbs_initial <- model.init$params$startProbs
#'         emissionParams_initial <- model.init$params$emissionParams
#'     }
#'   
#'     ### Define parameters for C function ###
#'     params <- list()
#'     params$startProbs_initial <- startProbs_initial
#'     params$transProbs_initial <- transProbs_initial
#'     params$emissionParams_initial <- emissionParams_initial
#'     params$transDist <- transDist
#'     params$eps <- eps
#'     params$maxtime <- max.time
#'     params$maxiter <- max.iter
#'     params$verbosity <- verbosity
#'     params$numThreads <- num.threads
#'     
#'     ### Call the C function ###
#'     on.exit(cleanup())
#'     if (is.null(fit.on.chrom)) {
#'         ptm <- startTimedMessage("Baum-Welch: Fitting HMM parameters\n")
#'         hmm <- fitHMM(counts=counts, distances=distances, params=params, algorithm=1)
#'         message("Time spent in Baum-Welch:", appendLF=FALSE)
#'         stopTimedMessage(ptm)
#'     } else {
#'         ptm <- startTimedMessage("Baum-Welch: Fitting HMM parameters\n")
#'         message(" ... on chromosomes ", paste0(fit.on.chrom, collapse=', '))
#'         hmm <- fitHMM(counts=counts, distances=distances, params=params, algorithm=1)
#'         message("Time spent in Baum-Welch:", appendLF=FALSE)
#'         stopTimedMessage(ptm)
#'         counts <- as.integer(data$observable)
#'         params2 <- list()
#'         params2$startProbs_initial <- hmm$startProbs
#'         params2$transProbs_initial <- hmm$transProbs
#'         params2$emissionParams_initial <- hmm$emissionParams
#'         params2$transDist <- transDist
#'         params2$eps <- eps
#'         params2$maxtime <- max.time
#'         params2$maxiter <- max.iter
#'         params2$verbosity <- verbosity
#'         params2$numThreads <- num.threads
#'         ptm <- startTimedMessage("Forward-Backward: Obtaining state sequence - no updates\n")
#'         message(" ... on all chromosomes")
#'         hmm <- fitHMM(counts=counts, distances=distances, params=params2, algorithm=2)
#'         message("Time spent in Forward-Backward:", appendLF=FALSE)
#'         stopTimedMessage(ptm)
#'     }
#'     
#'     ### Construct result object ###
#'     ptm <- startTimedMessage("Compiling results ...")
#'     r <- list()
#'     class(r) <- "NcomponentHMM"
#'     if (hmm$error == "") {
#'       
#'         r$convergenceInfo <- hmm$convergenceInfo
#'         names(hmm$weights) <- states
#'         r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, transDist=hmm$transDist, emissionParams=hmm$emissionParams, weights=hmm$weights)
#'         r$params.initial <- params
#'         # States and posteriors
#'         data$posteriorMax <- NA
#'         for (i1 in 0:(nrow(hmm$posteriors)-1)) {
#'             mask <- hmm$states == i1
#'             data$posteriorMax[mask] <- hmm$posteriors[i1+1,mask]
#'         }
#'         rownames(hmm$posteriors) <- states
#'         data$posteriors <- t(hmm$posteriors)
#'         data$state <- factor(states, levels=states)[hmm$states+1]
#'         ## Segmentation
#'         data.collapse <- data
#'         mcols(data.collapse) <- NULL
#'         data.collapse$state <- data$state
#'         df <- suppressMessages( collapseBins(as.data.frame(data.collapse), column2collapseBy = 'state') )
#'         segments <- methods::as(df, 'GRanges')
#'         seqlengths(segments) <- seqlengths(data)[seqlevels(segments)]
#'         
#'         ## Reorder the states by mean of the emission distrbution
#'         stateorder <- as.character(states[order(hmm$emissionParams$mu)])
#'         # Initial params
#'         r$params.initial$startProbs_initial <- r$params.initial$startProbs_initial[stateorder]
#'         names(r$params.initial$startProbs_initial) <- states
#'         r$params.initial$transProbs_initial <- r$params.initial$transProbs_initial[stateorder, stateorder]
#'         rownames(r$params.initial$transProbs_initial) <- states
#'         colnames(r$params.initial$transProbs_initial) <- states
#'         r$params.initial$emissionParams_initial <- r$params.initial$emissionParams_initial[stateorder, ]
#'         rownames(r$params.initial$emissionParams_initial) <- states
#'         # Params
#'         r$params$startProbs <- r$params$startProbs[stateorder]
#'         names(r$params$startProbs) <- states
#'         r$params$transProbs <- r$params$transProbs[stateorder, stateorder]
#'         rownames(r$params$transProbs) <- states
#'         colnames(r$params$transProbs) <- states
#'         r$params$emissionParams <- r$params$emissionParams[stateorder, ]
#'         rownames(r$params$emissionParams) <- states
#'         r$params$weights <- r$params$weights[stateorder]
#'         names(r$params$weights) <- states
#'         # Posteriors
#'         data$posteriors <- data$posteriors[,stateorder]
#'         colnames(data$posteriors) <- states
#'         # Actual states
#'         stateorder.mapping <- states
#'         names(stateorder.mapping) <- stateorder
#'         stateorder.mapping <- factor(stateorder.mapping, levels=states)
#'         data$state <- stateorder.mapping[as.character(data$state)]
#'         segments$state <- stateorder.mapping[as.character(segments$state)]
#'     }
#'     r$data <- data
#'     r$segments <- segments
#'     stopTimedMessage(ptm)
#'     
#'     return(r)
#' }
#' 
#' 
#' #' @importFrom stats pnbinom qnorm dnbinom
#' prepareMultivariate = function(data, hmms, statedef, max.states=Inf) {
#'   
#'     ## Assign variables
#'     bins <- data
#'     mcols(bins) <- NULL
#'     bins$distance <- data$distance
#'     modelnames <- names(hmms)
#'     l <- list()
#'     for (i1 in 1:ncol(statedef)) {
#'         l[[i1]] <- statedef[,i1]
#'     }
#'     states <- rownames(statedef)
#'     names(states) <- do.call(paste, l)
#'     mapping <- names(states)
#'     names(mapping) <- states
#'     components <- lapply(statedef, levels)
#'     if (any(names(components) != modelnames)) {
#'         stop("names(hmms) must be equal to names(statedef)")
#'     }
#'     numstates <- length(states)
#'     
#'     ## Go through HMMs and extract stuff
#'     ptm <- startTimedMessage("Extracting states ...")
#'     counts <- matrix(NA, nrow=length(bins), ncol=length(modelnames), dimnames=list(NULL, modelnames))
#'     emissionParams <- list()
#'     statevec <- list()
#'     for (i1 in 1:length(modelnames)) {
#'         modelname <- modelnames[i1]
#'         counts[,i1] <- hmms[[modelname]]$data$observable
#'         emissionParams[[modelname]] <- hmms[[modelname]]$params$emissionParams
#'         statevec[[modelname]] <- hmms[[modelname]]$data$state
#'     }
#'     maxcounts <- max(apply(counts, 2, max))
#'     bins$counts <- counts
#'     bins$state <- factor(states[do.call(paste, statevec)], levels=states)
#'     stopTimedMessage(ptm)
#'     
#'     # Clean up to reduce memory usage
#'     # remove(hmm)
#'     
#'     ## We pre-compute the z-values for each number of counts
#'     ptm <- startTimedMessage("Computing pre z-matrix ...")
#'     z.per.read <- list()
#'     for (modelname in modelnames) {
#'         z.per.read[[modelname]] <- array(NA, dim=c(maxcounts+1, length(components[[modelname]])), dimnames=list(counts=NULL, component=components[[modelname]]))
#'         xcounts = 0:maxcounts
#'         for (component in components[[modelname]]) {
#'             
#'             size = emissionParams[[modelname]][component,'size']
#'             prob = emissionParams[[modelname]][component,'prob']
#'             u = stats::pnbinom(xcounts, size, prob)
#'             
#'             # Check for infinities in u and set them to max value which is not infinity
#'             qnorm_u = stats::qnorm(u)
#'             mask <- qnorm_u==Inf
#'             qnorm_u[mask] <- stats::qnorm(1-1e-16)
#'             z.per.read[[modelname]][ , component] <- qnorm_u
#'           
#'         }
#'     }
#'     stopTimedMessage(ptm)
#'     
#'     ## Compute the z matrix
#'     ptm <- startTimedMessage("Transfering values into z-matrix ...")
#'     z.per.bin <- list()
#'     for (modelname in modelnames) {
#'         z.per.bin[[modelname]] = array(NA, dim=c(length(bins), length(components[[modelname]])), dimnames=list(bin=NULL, component=components[[modelname]]))
#'         for (component in components[[modelname]]) {
#'             z.per.bin[[modelname]][ , component] <- z.per.read[[modelname]][bins$counts[, modelname]+1, i1]
#'         }
#'     }
#'     
#'     # Clean up to reduce memory usage
#'     remove(z.per.read)
#'     stopTimedMessage(ptm)
#'     
#'     ### Calculate correlation matrix
#'     ptm <- startTimedMessage("Computing inverse of correlation matrix ...")
#'     correlationMatrix = array(NA, dim=c(length(modelnames),length(modelnames),length(states)), dimnames=list(model=modelnames, model=modelnames, state=states))
#'     for (i1 in 1:dim(correlationMatrix)[3]) {
#'         correlationMatrix[,,i1] <- diag(length(modelnames))
#'     }
#'     correlationMatrixInverse = correlationMatrix
#'     determinant = rep(1, length(states))
#'     names(determinant) <- states
#'     
#'     for (state in states) {
#'         istate = which(state==states)
#'         mask = which(bins$state==state)
#'         if (length(mask) > 100) {
#'             # Subselect z
#'             z.temp <- matrix(NA, nrow=length(mask), ncol=length(modelnames), dimnames=list(NULL, model=modelnames))
#'             for (modelname in modelnames) {
#'                 z.temp[,modelname] <- z.per.bin[[modelname]][mask, as.character(statedef[state, modelname])]
#'             }
#'             temp <- tryCatch({
#'                 correlationMatrix[,,state] <- cor(z.temp)
#'                 determinant[state] <- det( correlationMatrix[,,state] )
#'                 correlationMatrixInverse[,,state] <- chol2inv(chol(correlationMatrix[,,state]))
#'                 0
#'             }, warning = function(war) {
#'                 1
#'             }, error = function(err) {
#'                 1
#'             })
#'         }
#'         if (any(is.na(correlationMatrixInverse[,,state]))) {
#'             correlationMatrixInverse[,,state] <- diag(length(modelnames))
#'         }
#'     }
#'     remove(z.per.bin)
#'     stopTimedMessage(ptm)
#'     
#'     ### Put correlation matrix into list
#'     corrList <- list()
#'     corrInvList <- list()
#'     for (state in states) {
#'         corrList[[state]] <- correlationMatrix[,, state]
#'         corrInvList[[state]] <- correlationMatrixInverse[,, state]
#'     }
#'     correlationMatrix <- corrList
#'     correlationMatrixInverse <- corrInvList
#'     
#'     ### Select states to use ###
#'     weights <- sort(table(bins$state), decreasing = TRUE)
#'     states2use <- names(weights)[1:min(max.states, length(weights))]
#'     correlationMatrix <- correlationMatrix[states2use]
#'     correlationMatrixInverse <- correlationMatrixInverse[states2use]
#'     determinant <- determinant[states2use]
#'     statedef <- statedef[states2use,]
#'     mapping <- mapping[states2use]
#'     
#'     # Return parameters
#'     out <- list(
#'                 correlationMatrix = correlationMatrix,
#'                 correlationMatrixInverse = correlationMatrixInverse,
#'                 determinant = determinant,
#'                 observable = counts,
#'                 emissionParams = emissionParams,
#'                 statedef = statedef,
#'                 mapping = mapping,
#'                 weights = weights
#'     )
#'     return(out)
#' }


