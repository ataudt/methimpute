addDistance <- function(data) {
    ptm <- startTimedMessage("Adding distance ...")
    ## Define distance as "from previous to current"
    distance <- c(-1, start(data)[-1] - end(data)[-length(data)] - 1)
    distance[distance < 0] <- Inf 
    stopTimedMessage(ptm)
    return(distance)
}

addTransitionContext <- function(data) {
    ptm <- startTimedMessage("Adding transition context ...")
    
    ## Generating possible transition context levels
    contexts <- intersect(levels(data$context), unique(data$context))
    ncontexts <- length(contexts)
    transitionContexts <- character()
    for (c1 in 1:length(contexts)) {
        # for (c2 in c1:length(contexts)) {
        for (c2 in 1:length(contexts)) {
            transitionContexts[length(transitionContexts)+1] <- paste0(contexts[c1], '-', contexts[c2])
        }
    }
    ## Define transition context as "from previous to current"
    transitionContext <- factor(c(NA, paste0(data$context[-length(data)], '-', data$context[-1])), levels=transitionContexts)
    stopTimedMessage(ptm)
    return(transitionContext)
}