addDistance <- function(data, separate.contexts=FALSE) {
    ptm <- startTimedMessage("Adding distance ...")
    contexts <- intersect(levels(data$context), unique(data$context))
    ## Define distance as "from previous to current"
    if (separate.contexts) {
        distance <- numeric(length(data))
        for (context in contexts) {
            context.mask <- data$context == context
            distance[context.mask] <- c(-1, start(data)[context.mask][-1] - end(data)[context.mask][-length(which(context.mask))] - 1)
            distance[context.mask][distance[context.mask] < 0] <- Inf 
        }
    } else {
        distance <- c(-1, start(data)[-1] - end(data)[-length(data)] - 1)
        distance[distance < 0] <- Inf 
    }
    stopTimedMessage(ptm)
    return(distance)
}

addTransitionContext <- function(data, separate.contexts=FALSE) {
    ptm <- startTimedMessage("Adding transition context ...")
    
    ## Generating possible transition context levels
    contexts <- intersect(levels(data$context), unique(data$context))
    ncontexts <- length(contexts)
    transitionContexts <- character()
    for (c1 in 1:length(contexts)) {
        # for (c2 in c1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            if (separate.contexts & c1!=c2) {
                next
            }
            transitionContexts[length(transitionContexts)+1] <- paste0(contexts[c1], '-', contexts[c2])
        }
    }
    ## Define transition context as "from previous to current"
    if (separate.contexts) {
        transitionContext <- factor(NA, levels=transitionContexts)
        for (context in contexts) {
            context.mask <- data$context == context
            transitionContext[context.mask] <- factor(c(NA, paste0(data$context[context.mask][-length(which(context.mask))], '-', data$context[context.mask][-1])), levels=transitionContexts)
        }
    } else {
        transitionContext <- factor(c(NA, paste0(data$context[-length(data)], '-', data$context[-1])), levels=transitionContexts)
    }
    stopTimedMessage(ptm)
    return(transitionContext)
}