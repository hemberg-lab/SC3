#' @export
run_sc3.SCESet <- function(
                        object,
                        ks, 
                        exprs_values = "counts",
                        gene.filter = TRUE,
                        gene.filter.fraction = 0.06,
                        gene.reads.rare = 2,
                        gene.reads.ubiq = 0,
                        log.scale = TRUE,
                        d.region.min = 0.04,
                        d.region.max = 0.07,
                        k.means.iter.max = 1e+09,
                        k.means.nstart = 1000,
                        show.original.labels = FALSE,
                        svm.num.cells = NULL,
                        svm.train.inds = NULL,
                        n.cores = NULL,
                        seed = 1) {
    dataset <- object@assayData[[exprs_values]]
    if ( is.null(dataset) )
        warning(paste0("The object does not contain ", exprs_values, " expression values. Returning NULL."))
    res <- sc3(dataset = dataset,
        ks = ks,
        gene.filter = gene.filter,
        gene.filter.fraction = gene.filter.fraction,
        gene.reads.rare = gene.reads.rare,
        gene.reads.ubiq = gene.reads.ubiq,
        log.scale = log.scale,
        d.region.min = d.region.min,
        d.region.max = d.region.max,
        k.means.iter.max = k.means.iter.max,
        k.means.nstart = k.means.nstart,
        show.original.labels = show.original.labels,
        svm.num.cells = svm.num.cells,
        svm.train.inds = svm.train.inds,
        n.cores = n.cores,
        seed = seed)
    object@consensus <- res
    return(object)
}

#' @export
setMethod("run_sc3", signature(object = "SCESet"),
          function(
                object,
                ks, 
                exprs_values = "counts",
                gene.filter = TRUE,
                gene.filter.fraction = 0.06,
                gene.reads.rare = 2,
                gene.reads.ubiq = 0,
                log.scale = TRUE,
                d.region.min = 0.04,
                d.region.max = 0.07,
                k.means.iter.max = 1e+09,
                k.means.nstart = 1000,
                show.original.labels = FALSE,
                svm.num.cells = NULL,
                svm.train.inds = NULL,
                n.cores = NULL,
                seed = 1
          ) {
              run_sc3.SCESet(object,
                 ks,
                 exprs_values,
                 gene.filter,
                 gene.filter.fraction,
                 gene.reads.rare,
                 gene.reads.ubiq,
                 log.scale,
                 d.region.min,
                 d.region.max,
                 k.means.iter.max,
                 k.means.nstart,
                 show.original.labels,
                 svm.num.cells,
                 svm.train.inds,
                 n.cores,
                 seed)
          })


#' @export
estimate_k.SCESet <- function(
    object,
    exprs_values = "counts",
    gene.filter = TRUE,
    gene.filter.fraction = 0.06,
    gene.reads.rare = 2,
    gene.reads.ubiq = 0,
    log.scale = TRUE) {
    dataset <- object@assayData[[exprs_values]]
    if ( is.null(dataset) )
        warning(paste0("The object does not contain ", exprs_values, " expression values. Returning NULL."))
    res <- estkTW(dataset = dataset,
                  gene.filter = gene.filter,
                  gene.filter.fraction = gene.filter.fraction,
                  gene.reads.rare = gene.reads.rare,
                  gene.reads.ubiq = gene.reads.ubiq,
                  log.scale = log.scale)
    return(object)
}

#' @export
setMethod("estimate_k", signature(object = "SCESet"),
          function(
              object,
              exprs_values = "counts",
              gene.filter = TRUE,
              gene.filter.fraction = 0.06,
              gene.reads.rare = 2,
              gene.reads.ubiq = 0,
              log.scale = TRUE
          ) {
              estimate_k.SCESet(object,
                                exprs_values,
                                gene.filter,
                                gene.filter.fraction,
                                gene.reads.rare,
                                gene.reads.ubiq,
                                log.scale)
          })

