#' Generate a colour palette by k-means clustering of LAB colour space.
#'
#' Generate a palette of distinct colours through k-means clustering of LAB
#' colour space.
#'
#' @param n Numeric. The number of colours to generate.
#' @param hmin Numeric, in the range [0, 360]. The lower limit of the hue range
#'   to be clustered.
#' @param hmax Numeric, in the range [0, 360]. The upper limit of the hue range
#'   to be clustered.
#' @param cmin Numeric, in the range [0, 180]. The lower limit of the chroma
#'   range to be clustered.
#' @param cmax Numeric, in the range [0, 180]. The upper limit of the chroma
#'   range to be clustered.
#' @param lmin Numeric, in the range [0, 100]. The lower limit of the luminance
#'   range to be clustered.
#' @param lmax Numeric, in the range [0, 100]. The upper limit of the luminance
#'   range to be clustered.
#' @param plot Logical. Should the colour swatches be plotted (using
#'   \code{\link{swatch}})?
#' @param random Logical. If \code{TRUE}, clustering will be determined by the
#'   existing RNG state. If \code{FALSE}, the seed will be set to \code{1} for
#'   clustering, and on exit, the function will restore the pre-existing RNG
#'   state.
#' @return A vector of \code{n} colours (as hexadecimal strings), representing
#'   centers of clusters determined through k-means clustering of the LAB colour
#'   space delimited by \code{hmin}, \code{hmax}, \code{cmin}, \code{cmax},
#'   \code{lmin} and \code{lmax}.
#' @details Note that \code{iwanthue} currently doesn't support \code{hmin}
#'   greater than \code{hmax} (which should be allowed, since hue is circular).
#' @references
#' \itemize{
#'   \item \href{https://github.com/johnbaums/hues}{R implementation of
#'   iwanthue by John Baumgartner}
#'   \item \href{http://tools.medialab.sciences-po.fr/iwanthue/}{iwanthue -
#'   colors for data scientists}
#'   \item \href{https://github.com/medialab/iwanthue}{iwanthue on
#'   GitHub}
#' }
#' @seealso \code{\link{swatch}}
#' @importFrom colorspace LAB hex coords
#' @importFrom methods as
#' @importFrom stats kmeans
iwanthue <- function(n, hmin = 0, hmax = 360, cmin = 0, cmax = 180, lmin = 0, lmax = 100, 
    plot = FALSE, random = FALSE) {
    stopifnot(hmin >= 0, cmin >= 0, lmin >= 0, hmax <= 360, cmax <= 180, lmax <= 
        100, hmin <= hmax, cmin <= cmax, lmin <= lmax, n > 0)
    if (!random) {
        if (exists(".Random.seed", .GlobalEnv)) {
            old_seed <- .GlobalEnv$.Random.seed
            on.exit(.GlobalEnv$.Random.seed <- old_seed)
        } else {
            on.exit(rm(".Random.seed", envir = .GlobalEnv))
        }
        set.seed(1)
    }
    lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1), seq(-100, 100, 5), seq(-110, 
        100, 5))))
    if (any((hmin != 0 || cmin != 0 || lmin != 0 || hmax != 360 || cmax != 180 || 
        lmax != 100))) {
        hcl <- methods::as(lab, "polarLUV")
        hcl_coords <- coords(hcl)
        hcl <- hcl[which(hcl_coords[, "H"] <= hmax & hcl_coords[, "H"] >= hmin & 
            hcl_coords[, "C"] <= cmax & hcl_coords[, "C"] >= cmin & hcl_coords[, 
            "L"] <= lmax & hcl_coords[, "L"] >= lmin), ]
        lab <- methods::as(hcl, "LAB")
    }
    lab <- lab[which(!is.na(hex(lab))), ]
    clus <- kmeans(coords(lab), n, iter.max = 50)
    if (isTRUE(plot)) {
        swatch(hex(LAB(clus$centers)))
    }
    hex(LAB(clus$centers))
}

#' Plot colour swatches for a vector of colours
#'
#' Plot named colour swatches for a vector of colours.
#'
#' @param x a vector of colours, specified as: colour names (i.e.
#' colour names returned by \code{colors()}); numeric indices into
#' \code{palette()}, or hexadecimal strings in the form \code{'#RRGGBB'}, where
#' \code{RR}, \code{GG}, and \code{BB} are pairs of hexadecimal digits
#' representing red, green, and blue components, in the range \code{00} to
#' \code{FF}.
#' @return \code{NULL}. The colour swatch is plotted to the active plotting
#'   device.
#' @seealso \code{\link{iwanthue}}
#'
#' @importFrom graphics par strwidth barplot
#'
swatch <- function(x) {
    par(mai = c(0.2, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.2, 0.4))
    barplot(rep(1, length(x)), col = rev(x), space = 0.1, axes = FALSE, names.arg = rev(x), 
        cex.names = 0.8, horiz = TRUE, las = 1)
    return(invisible(NULL))
}
