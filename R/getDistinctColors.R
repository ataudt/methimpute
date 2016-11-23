

#' Get distinct colors
#' 
#' Get a set of distinct colors selected from \code{\link{colors}}.
#' 
#' The function computes the euclidian distance between all \code{\link{colors}} and iteratively selects those that have the furthest closes distance to the set of already selected colors.
#' 
#' @param n Number of colors to select. If \code{n} is a character vector, \code{length(n)} will be taken as the number of colors and the colors will be named by \code{n}.
#' @param start.color Color to start the selection process from.
#' @param exclude.colors Character vector with colors that should not be used.
#' @param exclude.rgb.above Exclude colors where all RGB values are above. This is useful to exclude whitish colors.
#' @return A character vector with colors.
#' @author Aaron Taudt
#' @importFrom grDevices col2rgb
#' @importFrom stats dist
#' @export
#' @examples
#'cols <- getDistinctColors(5)
#'pie(rep(1,5), labels=cols, col=cols)
#'
getDistinctColors <- function(n, start.color='blue4', exclude.colors=c('white','black','gray','grey'), exclude.rgb.above=210) {
    
    if (is.character(n)) {
        n.names <- n
        n <- length(n)
    }
    cols <- grDevices::colors()
    # Exclude unwanted colors
    cols <- grep(paste(exclude.colors, collapse='|'), cols, invert=TRUE, value=TRUE)
    # Get RGB values
    rgbs <- t(grDevices::col2rgb(cols))
    rownames(rgbs) <- cols
    # Exclude whitish colors
    rgbs <- rgbs[apply(rgbs, 1, function(x) { !all(x>exclude.rgb.above) }), ]
    # Calculate distance
    coldist <- as.matrix(stats::dist(rgbs, method='euclidean'))
    
    if (n == 1) {
        return(start.color)
    }
    # Iteratively select colors
    selected.cols <- character()
    selected.cols[1] <- start.color
    for (i1 in 2:n) {
        m <- as.matrix(coldist[,selected.cols])
        closest.dist <- apply(m, 1, min)
        furthest.dist <- which.max(closest.dist)
        selected.cols[i1] <- names(furthest.dist)
        if (selected.cols[i1] == selected.cols[i1-1]) {
            selected.cols <- selected.cols[-length(selected.cols)]
            break
        }
    }
    colors <- rep(selected.cols, ceiling(n / length(selected.cols)))[1:n]
    if (length(colors) > length(selected.cols)) {
        warning("Recycling colors because we only have ", length(selected.cols), " distinct colors to choose from.")
    }
    names(colors) <- n.names
    return(colors)
    
}
