#' @export
setGeneric("sc3", signature = "object", function(object, ks = NULL, 
        gene_filter = TRUE, pct_dropout_min = 10, pct_dropout_max = 90, 
        d_region_min = 0.04, d_region_max = 0.07, svm_num_cells = NULL, 
        svm_train_inds = NULL, svm_max = 5000, n_cores = NULL, kmeans_nstart = NULL, 
        kmeans_iter_max = 1e+09, k_estimator = FALSE, biology = FALSE, rand_seed = 1) {
    standardGeneric("sc3")
})

#' @export
setGeneric("sc3_estimate_k", signature = "object", function(object) {
    standardGeneric("sc3_estimate_k")
})

#' @export
setGeneric("sc3_prepare", function(object, gene_filter = TRUE, 
        pct_dropout_min = 10, pct_dropout_max = 90, d_region_min = 0.04, 
        d_region_max = 0.07, svm_num_cells = NULL, svm_train_inds = NULL, 
        svm_max = 5000, n_cores = NULL, kmeans_nstart = NULL, 
        kmeans_iter_max = 1e+09, rand_seed = 1) {
    standardGeneric("sc3_prepare")
})

#' @export
setGeneric("sc3_calc_dists", signature = "object", function(object) {
    standardGeneric("sc3_calc_dists")
})

#' @export
setGeneric("sc3_calc_transfs", signature = "object", function(object) {
    standardGeneric("sc3_calc_transfs")
})

#' @export
setGeneric("sc3_kmeans", signature = "object", function(object, ks = NULL) {
    standardGeneric("sc3_kmeans")
})

#' @export
setGeneric("sc3_calc_consens", signature = "object", function(object) {
    standardGeneric("sc3_calc_consens")
})

#' @export
setGeneric("sc3_calc_biology", signature = "object", function(object, ks = NULL, 
                                                              regime = NULL) {
    standardGeneric("sc3_calc_biology")
})

#' @export
setGeneric("sc3_interactive", signature = "object", function(object) {
    standardGeneric("sc3_interactive")
})

#' @export
setGeneric("sc3_run_svm", signature = "object", function(object, ks = NULL) {
    standardGeneric("sc3_run_svm")
})

#' @export
setGeneric("sc3_plot_consensus", signature = "object", function(object, k, 
                                                                show_pdata = NULL) {
    standardGeneric("sc3_plot_consensus")
})

#' @export
setGeneric("sc3_plot_silhouette", signature = "object", function(object, k) {
    standardGeneric("sc3_plot_silhouette")
})

#' @export
setGeneric("sc3_plot_expression", signature = "object", function(object, k, show_pdata = NULL) {
    standardGeneric("sc3_plot_expression")
})

#' @export
setGeneric("sc3_plot_de_genes", signature = "object", function(object, 
                                        k, p.val = 0.01, show_pdata = NULL) {
    standardGeneric("sc3_plot_de_genes")
})

#' @export
setGeneric("sc3_plot_markers", signature = "object", function(object, k, auroc = 0.85, 
                                        p.val = 0.01, show_pdata = NULL) {
    standardGeneric("sc3_plot_markers")
})

#' @export
setGeneric("sc3_plot_cluster_stability", signature = "object", function(object, k) {
    standardGeneric("sc3_plot_cluster_stability")
})

#' @export
setGeneric("sc3_export_results_xls", signature = "object", function(object, 
                                                filename = "sc3_results.xls") {
    standardGeneric("sc3_export_results_xls")
})
