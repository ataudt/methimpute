#' Call methylated regions
#' 
#' Call methylated, unmethylated and hemimethylated regions by fitting a multivariate Hidden Markov Model to the data.
#' 
#' In a first step, the function will use \code{\link{fitSignalBackground}} to obtain the distribution parameters for a signal-background model for methylated and unmethylated read counts, respectively. The resulting distributions are used to fit a multivariate HMM with states \code{c("Background", "Methylated", "UNmethylated", "Hemimethylated")}.
#' 
#' @inheritParams fitSignalBackground
#' @return A list with fitted parameters, posteriors and input parameters.
#' 
callMethylation <- function(data, fit.on.chrom=NULL, transDist=10000, eps=0.01, max.time=Inf, max.iter=Inf, count.cutoff=1000, verbosity=1) {
  
    ### Input checks ###
    if (is.null(max.time)) {
        max.time <- Inf
    }
    if (is.null(max.iter)) {
        max.iter <- Inf
    }
  
    ### Assign variables ###
    modelnames <- c('minus', 'plus')
    components <- c('background', 'signal')
    states <- c("Background", "Methylated", "UNmethylated", "Hemimethylated")
    numstates <- length(states)
    
    ### State definition and mapping ###
    statedef <- permute(components, modelnames)
    rownames(statedef) <- states
    
    ### Fit HMMs for p- and m-counts ###
    messageU("Fitting two-component HMM", underline="=", overline='=')
    mhmm <- fitSignalBackground(data, observable='mcounts', fit.on.chrom=fit.on.chrom, transDist=transDist, eps=eps, max.time=max.time, max.iter=max.iter, count.cutoff=count.cutoff, verbosity=verbosity)
    phmm <- fitSignalBackground(data, observable='pcounts', fit.on.chrom=fit.on.chrom, transDist=transDist, eps=eps, max.time=max.time, max.iter=max.iter, count.cutoff=count.cutoff, verbosity=verbosity)
    hmms <- list(mhmm, phmm)
    names(hmms) <- modelnames
    
    ### Initial probabilities ###
    s <- 0.9
    transProbs_initial <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=states, to=states))
    diag(transProbs_initial) <- s
    startProbs_initial <- rep(1/numstates, numstates)
    names(startProbs_initial) <- states
    
    ### Prepare the bivariate HMM ###
    messageU("Multivariate HMM", underline="=", overline="=")
    cormat <- prepareMultivariate(data=data, hmms=hmms, statedef=statedef)
    data$observable <- cormat$observable
    
    ### Define parameters for C function ###
    params <- list()
    params$startProbs_initial <- startProbs_initial
    params$transProbs_initial <- transProbs_initial
    params$emissionParamsList <- cormat$emissionParams
    params$transDist <- transDist
    params$eps <- eps
    params$maxtime <- max.time
    params$maxiter <- max.iter
    params$verbosity <- verbosity
    params$correlationMatrixInverse <- cormat$correlationMatrixInverse
    params$determinant <- cormat$determinant
    params$statedef <- statedef
    
    ### Fit the HMM ###
    message("Baum-Welch: Fitting parameters for multivariate HMM")
    on.exit(cleanup())
    hmm <- fitMultiHMM(data$observable, data$distance, params)
    
    ### Construct result object ###
    ptm <- startTimedMessage("Compiling results ...")
    r <- list()
    if (hmm$error == "") {
        r$convergenceInfo <- hmm$convergenceInfo
        names(hmm$weights) <- states
        r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, emissionParamsList=params$emissionParamsList, weights=hmm$weights)
        r$params.initial <- params
        # States and posteriors
        rownames(hmm$posteriors) <- states
        data$posteriors <- t(hmm$posteriors)
        data$state <- factor(states, levels=states)[hmm$states+1]
    }
    r$data <- data
    stopTimedMessage(ptm)
    
    return(r)
}



#' Fit a two-component HMM
#' 
#' Fit a two-component Hidden Markov Model to the supplied counts. The transition matrix is distance-dependent with exponential decaying constant \code{transDist}. Components are modeled as negative binomial distributions.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'distance' and 'counts'.
#' @param observable A character naming the column of \code{data} that will serve as observable for the HMM.
#' @param fit.on.chrom A character vector giving the chromosomes on which the HMM will be fitted.
#' @param transDist The exponential decaying constant for the distance-dependent transition matrix. Should be given in the same units as \code{distances}.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param max.time Maximum running time in seconds for the Baum-Welch algorithm.
#' @param max.iter Maximum number of iterations for the Baum-Welch algorithm.
#' @param count.cutoff \code{counts} > \code{count.cutoff} will be set to \code{count.cutoff}. The purpose is to increase the speed of the Baum-Welch algorithm, which is slowed down by huge values in \code{counts}.
#' @param verbosity Integer from c(0,1) specifying the verbosity of the fitting procedure.
#' @return A list with fitted parameters, posteriors, and the input parameters.
#' 
fitSignalBackground <- function(data, observable='counts', fit.on.chrom=NULL, transDist=10000, eps=0.01, max.time=Inf, max.iter=Inf, count.cutoff=1000, verbosity=1) {
  
    ### Input checks ###
    if (is.null(max.time)) {
        max.time <- Inf
    }
    if (is.null(max.iter)) {
        max.iter <- Inf
    }
  
    ### Assign variables ###
    states <- c("background", "signal")
    numstates <- length(states)
    counts <- mcols(data)[,observable]
    distances <- data$distance
    
    ## Filter counts by cutoff
    counts[counts > count.cutoff] <- count.cutoff
    data$observable <- counts # assign it now to have filtered values in there
    
    ## Subset by chromosomes
    if (!is.null(fit.on.chrom)) {
        counts <- counts[as.logical(data@seqnames %in% fit.on.chrom)]
    }
  
    ### Initial probabilities ###
    s <- 0.9
    transProbs_initial <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=states, to=states))
    diag(transProbs_initial) <- s
    startProbs_initial <- rep(1/numstates, numstates)
    names(startProbs_initial) <- states
    
    ### Initialization of emission distributions ###
    counts.greater.0 <- counts[counts>0]
    mean.counts <- mean(counts.greater.0)
    var.counts <- var(counts.greater.0)
    ep <- list()
    ep[["background"]] <- data.frame(type='dnbinom', mu=mean.counts, var=var.counts)
    ep[["signal"]] <- data.frame(type='dnbinom', mu=mean.counts*1.5, var=var.counts*2)
    ep <- do.call(rbind, ep)
    # Make sure variance is bigger than mean for negative binomial
    ep$var[ep$mu >= ep$var] <- ep$mu + 1
    ep$size <- dnbinom.size(ep$mu, ep$var)
    ep$prob <- dnbinom.prob(ep$mu, ep$var)
    emissionParams_initial <- ep
  
    ### Define parameters for C function ###
    params <- list()
    params$startProbs_initial <- startProbs_initial
    params$transProbs_initial <- transProbs_initial
    params$emissionParams_initial <- emissionParams_initial
    params$transDist <- transDist
    params$eps <- eps
    params$maxtime <- max.time
    params$maxiter <- max.iter
    params$verbosity <- verbosity
    
    ### Call the C function ###
    on.exit(cleanup())
    if (is.null(fit.on.chrom)) {
        message("Baum-Welch: Fitting HMM parameters")
        hmm <- fitHMM(counts=counts, distances=distances, params=params, algorithm=1)
    } else {
        message("Baum-Welch: Fitting HMM parameters")
        message(" ... on chromosomes ", paste0(fit.on.chrom, collapse=', '))
        hmm <- fitHMM(counts=counts, distances=distances, params=params, algorithm=1)
        counts <- data$observable
        params2 <- list()
        params2$startProbs_initial <- hmm$startProbs
        params2$transProbs_initial <- hmm$transProbs
        params2$emissionParams_initial <- hmm$emissionParams
        params2$transDist <- transDist
        params2$eps <- eps
        params2$maxtime <- max.time
        params2$maxiter <- max.iter
        params2$verbosity <- verbosity
        message("Viterbi: Obtaining state sequence")
        message(" ... on all chromosomes")
        hmm <- fitHMM(counts=counts, distances=distances, params=params2, algorithm=2)
    }
    
    ### Construct result object ###
    ptm <- startTimedMessage("Compiling results ...")
    r <- list()
    if (hmm$error == "") {
        r$convergenceInfo <- hmm$convergenceInfo
        names(hmm$weights) <- states
        r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, emissionParams=hmm$emissionParams, weights=hmm$weights)
        r$params.initial <- params
        # States and posteriors
        rownames(hmm$posteriors) <- states
        data$posteriors <- t(hmm$posteriors)
        data$state <- factor(states, levels=states)[hmm$states+1]
    }
    r$data <- data
    stopTimedMessage(ptm)
    
    return(r)
}


#' Fit a three-component HMM for methylation levels
#' 
#' Fit a three-component Hidden Markov Model to the supplied counts. The transition matrix is distance-dependent with exponential decaying constant \code{transDist}. Components are modeled as beta distributions.
#' 
#' @param data A \code{\link[GenomicRanges]{GRanges}} object with metadata columns 'distance' and 'ratio'.
#' @param fit.on.chrom A character vector giving the chromosomes on which the HMM will be fitted.
#' @param transDist The exponential decaying constant for the distance-dependent transition matrix. Should be given in the same units as \code{distances}.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param max.time Maximum running time in seconds for the Baum-Welch algorithm.
#' @param max.iter Maximum number of iterations for the Baum-Welch algorithm.
#' @param verbosity Integer from c(0,1) specifying the verbosity of the fitting procedure.
#' @return A list with fitted parameters, posteriors, and the input parameters.
#' 
fitRatio <- function(data, fit.on.chrom=NULL, transDist=10000, eps=0.01, max.time=Inf, max.iter=Inf, verbosity=1) {
  
    ### Input checks ###
    if (is.null(max.time)) {
        max.time <- Inf
    }
    if (is.null(max.iter)) {
        max.iter <- Inf
    }
  
    ### Assign variables ###
    states <- c("UNmethylated", "Hemimethylated", "Methylated")
    numstates <- length(states)
    ratio <- data$ratio
    distances <- data$distance
    
    ## Subset by chromosomes
    if (!is.null(fit.on.chrom)) {
        ratio <- ratio[as.logical(data@seqnames %in% fit.on.chrom)]
    }
  
    ### Initial probabilities ###
    s <- 0.9
    transProbs_initial <- matrix((1-s)/(numstates-1), ncol=numstates, nrow=numstates, dimnames=list(from=states, to=states))
    diag(transProbs_initial) <- s
    startProbs_initial <- rep(1/numstates, numstates)
    names(startProbs_initial) <- states
    
    ### Initialization of emission distributions ###
    ep <- list()
    ep[["UNmethylated"]] <- data.frame(type='dbeta', a=1, b=4)
    ep[["Hemimethylated"]] <- data.frame(type='dbeta', a=10, b=10)
    ep[["Methylated"]] <- data.frame(type='dbeta', a=4, b=1)
    ep <- do.call(rbind, ep)
    emissionParams_initial <- ep
  
    ### Define parameters for C function ###
    params <- list()
    params$startProbs_initial <- startProbs_initial
    params$transProbs_initial <- transProbs_initial
    params$emissionParams_initial <- emissionParams_initial
    params$transDist <- transDist
    params$eps <- eps
    params$maxtime <- max.time
    params$maxiter <- max.iter
    params$verbosity <- verbosity
    
    ### Call the C function ###
    on.exit(cleanup())
    if (is.null(fit.on.chrom)) {
        message("Baum-Welch: Fitting HMM parameters")
        hmm <- fitHMMratio(ratio=ratio, distances=distances, params=params, algorithm=1)
    } else {
        message("Baum-Welch: Fitting HMM parameters")
        message(" ... on chromosomes ", paste0(fit.on.chrom, collapse=', '))
        hmm <- fitHMMratio(ratio=ratio, distances=distances, params=params, algorithm=1)
        ratio <- data$ratio
        params2 <- list()
        params2$startProbs_initial <- hmm$startProbs
        params2$transProbs_initial <- hmm$transProbs
        params2$emissionParams_initial <- hmm$emissionParams
        params2$transDist <- transDist
        params2$eps <- eps
        params2$maxtime <- max.time
        params2$maxiter <- max.iter
        params2$verbosity <- verbosity
        message("Viterbi: Obtaining state sequence")
        message(" ... on all chromosomes")
        hmm <- fitHMMratio(ratio=ratio, distances=distances, params=params2, algorithm=2)
    }
    
    ### Construct result object ###
    ptm <- startTimedMessage("Compiling results ...")
    r <- list()
    if (hmm$error == "") {
        r$convergenceInfo <- hmm$convergenceInfo
        names(hmm$weights) <- states
        r$params <- list(startProbs=hmm$startProbs, transProbs=hmm$transProbs, emissionParams=hmm$emissionParams, weights=hmm$weights)
        r$params.initial <- params
        # States and posteriors
        rownames(hmm$posteriors) <- states
        data$posteriors <- t(hmm$posteriors)
        rownames(hmm$densities) <- states
        data$densities <- t(hmm$densities)
        data$state <- factor(states, levels=states)[hmm$states+1]
    }
    r$data <- data
    stopTimedMessage(ptm)
    
    return(r)
}


#' @importFrom stats pnbinom qnorm dnbinom
prepareMultivariate = function(data, hmms, statedef) {
  
    ## Assign variables
    bins <- data
    mcols(bins) <- NULL
    bins$distance <- data$distance
    modelnames <- names(hmms)
    l <- list()
    for (i1 in 1:ncol(statedef)) {
        l[[i1]] <- statedef[,i1]
    }
    states <- rownames(statedef)
    names(states) <- do.call(paste, l)
    components <- levels(statedef[,1])
    numstates <- length(states)
    
    ## Go through HMMs and extract stuff
    ptm <- startTimedMessage("Extracting states ...")
    counts <- matrix(NA, nrow=length(bins), ncol=2, dimnames=list(NULL, modelnames))
    emissionParams <- list()
    statevec <- list()
    for (i1 in 1:length(modelnames)) {
        modelname <- modelnames[i1]
        counts[,i1] <- hmms[[modelname]]$data$observable
        emissionParams[[modelname]] <- hmms[[modelname]]$params$emissionParams
        statevec[[modelname]] <- hmms[[modelname]]$data$state
    }
    maxcounts <- max(apply(counts, 2, max))
    bins$counts <- counts
    bins$state <- factor(states[do.call(paste, statevec)], levels=states)
    stopTimedMessage(ptm)
    
    # Clean up to reduce memory usage
    # remove(hmm)
    
    ## We pre-compute the z-values for each number of counts
    ptm <- startTimedMessage("Computing pre z-matrix ...")
    z.per.read <- array(NA, dim=c(maxcounts+1, length(modelnames), length(components)), dimnames=list(counts=NULL, model=modelnames, component=components))
    xcounts = 0:maxcounts
    for (modelname in modelnames) {
        for (component in components) {
            
            size = emissionParams[[modelname]][component,'size']
            prob = emissionParams[[modelname]][component,'prob']
            u = stats::pnbinom(xcounts, size, prob)
            
            # Check for infinities in u and set them to max value which is not infinity
            qnorm_u = stats::qnorm(u)
            mask <- qnorm_u==Inf
            qnorm_u[mask] <- stats::qnorm(1-1e-16)
            z.per.read[ , modelname, component] <- qnorm_u
          
        }
    }
    stopTimedMessage(ptm)
    
    ## Compute the z matrix
    ptm <- startTimedMessage("Transfering values into z-matrix ...")
    z.per.bin = array(NA, dim=c(length(bins), length(modelnames), length(components)), dimnames=list(bin=NULL, model=modelnames, component=components))
    for (modelname in modelnames) {
        for (component in components) {
            z.per.bin[ , modelname, component] <- z.per.read[bins$counts[, modelname]+1, modelname, i1]
        }
    }
    
    # Clean up to reduce memory usage
    remove(z.per.read)
    stopTimedMessage(ptm)
    
    ### Calculate correlation matrix
    ptm <- startTimedMessage("Computing inverse of correlation matrix ...")
    correlationMatrix = array(NA, dim=c(length(modelnames),length(modelnames),length(states)), dimnames=list(model=modelnames, model=modelnames, state=states))
    correlationMatrixInverse = correlationMatrix
    determinant = rep(NA, length(states))
    names(determinant) <- states
    
    for (state in states) {
        istate = which(state==states)
        mask = which(bins$state==state)
        # Subselect z
        z.temp <- matrix(NA, nrow=length(mask), ncol=length(modelnames), dimnames=list(NULL, model=modelnames))
        for (modelname in modelnames) {
            z.temp[,modelname] <- z.per.bin[mask, modelname, statedef[state, modelname]]
        }
        temp <- tryCatch({
            if (nrow(z.temp) > 100) {
                correlationMatrix[,,state] <- cor(z.temp)
                determinant[state] <- det( correlationMatrix[,,state] )
                correlationMatrixInverse[,,state] <- chol2inv(chol(correlationMatrix[,,state]))
            } else {
                correlationMatrix[,,state] <- diag(length(modelnames))
                determinant[state] <- 1 # det(diag(x)) == 1
                correlationMatrixInverse[,,state] <- diag(length(modelnames)) # solve(diag(x)) == diag(x)
            }
            0
        }, warning = function(war) {
            1
        }, error = function(err) {
            1
        })
        if (temp!=0) {
            correlationMatrix[,,state] <- diag(length(modelnames))
            determinant[state] <- 1
            correlationMatrixInverse[,,state] <- diag(length(modelnames))
        }
        if (any(is.na(correlationMatrixInverse[,,state]))) {
            correlationMatrixInverse[,,state] <- diag(length(modelnames))
        }
    }
    remove(z.per.bin)
    stopTimedMessage(ptm)
    
    ### Put correlation matrix into list
    corrList <- list()
    corrInvList <- list()
    for (state in states) {
        corrList[[state]] <- correlationMatrix[,, state]
        corrInvList[[state]] <- correlationMatrixInverse[,, state]
    }
    correlationMatrix <- corrList
    correlationMatrixInverse <- corrInvList
    
    # Return parameters
    out <- list(
                correlationMatrix = correlationMatrix,
                correlationMatrixInverse = correlationMatrixInverse,
                determinant = determinant,
                observable = counts,
                emissionParams = emissionParams
    )
    return(out)
}


