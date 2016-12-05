permute <- function(components) {
    ncol <- length(components)
    ncomponents <- c(sapply(components, length), 1)
    total.rows <- prod(ncomponents)
    perms <- list()
    for (i1 in ncol:1) {
        p <- prod(ncomponents[(i1+1):length(ncomponents)])
        pp1 <- prod(ncomponents[i1:length(components)])
        perms[[i1]] <- rep(rep(components[[i1]], each = p), total.rows / pp1)
    }
    perms <- as.data.frame(perms)
    names(perms) <- names(components)
    return(perms)
}