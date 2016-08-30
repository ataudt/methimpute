permute <- function(components, modelnames) {
    ncomponents <- length(components)
    ncol <- length(modelnames)
    # perms <- matrix(NA, nrow=ncomponents^ncol, ncol=ncol)
    perms <- list()
    for (i1 in ncol:1) {
        perms[[i1]] <- rep(rep(components, each=ncomponents^(ncol-i1)), ncomponents^(i1-1))
    }
    perms <- as.data.frame(perms)
    names(perms) <- modelnames
    return(perms)
}