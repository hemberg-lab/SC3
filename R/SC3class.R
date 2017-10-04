#' SC3 class
#' 
#' SC3 class
#' 
#' @import SingleCellExperiment
#' @export
setClass("SC3class",
         slots = list(sc3 = "list"),
         contains = "SingleCellExperiment"
)

#' Create a new object of SC3 class
#' 
#' Create a new object of SC3 class
#' 
#' @param object an object of `SingleCellExperiment` class
#' 
#' @return an object of `SC3` class
#' 
#' @importFrom methods slot slot<- slotNames
#' @export
newSC3 <- function(object) {
    d <- new("SC3class")
    for(s in slotNames(object)) {
        slot(d, s) <- slot(object, s)
    }
    return(d)
}
