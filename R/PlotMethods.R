#' Plot consensus matrix as a heatmap
#' 
#' The consensus matrix is a NxN 
#' matrix, where N is the number of cells.
#' It represents similarity between the cells based 
#' on the averaging of clustering results from all 
#' combinations of clustering parameters. Similarity 0 
#' (blue) means that the two cells are always assigned to different clusters. 
#' In contrast, similarity 1 (red) means that the two cells are always assigned 
#' to the same cluster. The consensus matrix is clustered by hierarchical 
#' clustering and has a diagonal-block structure. Intuitively, the perfect 
#' clustering is achieved when all diagonal blocks are completely red 
#' and all off-diagonal elements are completely blue.
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
#' 
#' @export
sc3_plot_consensus.SCESet <- function(object, k, show_pdata = NULL) {
    
    res <- prepare_output(object, k)
    
    consensus <- object@sc3$consensus[[as.character(k)]]$consensus
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(consensus)
        }
    }
    do.call(pheatmap::pheatmap, c(list(consensus, cluster_rows = res$hc, cluster_cols = res$hc, 
        cutree_rows = k, cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), 
        list(annotation_col = ann)[add_ann_col]))
}

#' @rdname sc3_plot_consensus.SCESet
#' @aliases sc3_plot_consensus
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_plot_consensus", signature(object = "SCESet"), function(object, k, 
    show_pdata = NULL) {
    sc3_plot_consensus.SCESet(object, k, show_pdata)
})

#' Plot silhouette indexes of the cells
#' 
#' A silhouette is a quantitative measure of the diagonality of the consensus 
#' matrix. An average silhouette width (shown at the bottom left of the silhouette 
#' plot) varies from 0 to 1, where 1 represents a perfectly block-diagonal 
#' consensus matrix and 0 represents a situation where there is no 
#' block-diagonal structure. The best clustering is achieved when the average 
#' silhouette width is close to 1.
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' 
#' @export
sc3_plot_silhouette.SCESet <- function(object, k) {
    res <- prepare_output(object, k)
    plot(res$silh, col = "black")
}

#' @rdname sc3_plot_silhouette.SCESet
#' @aliases sc3_plot_silhouette
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_plot_silhouette", signature(object = "SCESet"), function(object, k) {
    sc3_plot_silhouette.SCESet(object, k)
})

#' Plot expression matrix used for SC3 clustering as a heatmap
#' 
#' The expression panel represents the original input expression matrix 
#' (cells in columns and genes in rows) after cell and gene filters. 
#' Genes are clustered by kmeans with k = 100 (dendrogram on the left) and 
#' the heatmap represents the expression levels of the gene cluster centers 
#' after log2-scaling.
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
#' 
#' @export
sc3_plot_expression.SCESet <- function(object, k, show_pdata = NULL) {
    
    res <- prepare_output(object, k)
    
    dataset <- get_processed_dataset(object)
    
    if (!is.null(object@sc3$svm_train_inds)) {
        dataset <- dataset[, object@sc3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    do.call(pheatmap::pheatmap, c(list(dataset, cluster_cols = res$hc, kmeans_k = 100, 
        cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname sc3_plot_expression.SCESet
#' @aliases sc3_plot_expression
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_plot_expression", signature(object = "SCESet"), function(object, k, 
    show_pdata = NULL) {
    sc3_plot_expression.SCESet(object, k, show_pdata)
})

#' Plot expression of DE genes of the clusters identified by \code{SC3} as a heatmap
#' 
#' \code{SC3} plots gene expression profiles of the 50 genes with the lowest p-values. 
#' 
#' @name sc3_plot_de_genes
#' @aliases sc3_plot_de_genes, sc3_plot_de_genes,SCESet-method
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
#' 
#' @export
sc3_plot_de_genes.SCESet <- function(object, k, p.val = 0.01, show_pdata = NULL) {
    res <- prepare_output(object, k)
    
    dataset <- get_processed_dataset(object)
    
    if (!is.null(object@sc3$svm_train_inds)) {
        dataset <- dataset[, object@sc3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    de_genes <- organise_de_genes(object, k, p.val)
    de_genes <- head(de_genes, 50)
    row_ann <- data.frame(log10_padj = -log10(de_genes))
    rownames(row_ann) <- names(de_genes)
    
    do.call(pheatmap::pheatmap, c(list(dataset[names(de_genes), , drop = FALSE], show_colnames = FALSE, 
        cluster_rows = FALSE, cluster_cols = res$hc, cutree_cols = k, annotation_row = row_ann, 
        cellheight = 10), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname sc3_plot_de_genes
#' @aliases sc3_plot_de_genes
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_plot_de_genes", signature(object = "SCESet"), function(object, k, 
    p.val = 0.01, show_pdata = NULL) {
    sc3_plot_de_genes.SCESet(object, k, p.val, show_pdata)
})


#' Plot expression of marker genes of the clusters identified by SC3 as a heatmap
#' 
#' To find marker genes, for each gene a binary classifier is constructed 
#' based on the mean cluster expression values. The classifier prediction 
#' is then calculated using the gene expression ranks. The area under the 
#' receiver operating characteristic (ROC) curve is used to quantify the accuracy 
#' of the prediction. A p-value is assigned to each gene by using the Wilcoxon 
#' signed rank test. By default the genes with the area under the ROC curve (AUROC) > 0.85 
#' and with the p-value < 0.01 are selected and the top 10 marker 
#' genes of each cluster are visualized in this heatmap.
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' @param auroc area under the ROC curve
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
#' 
#' @export
sc3_plot_markers.SCESet <- function(object, k, auroc = 0.85, p.val = 0.01, show_pdata = NULL) {
    
    res <- prepare_output(object, k)
    
    dataset <- get_processed_dataset(object)
    
    if (!is.null(object@sc3$svm_train_inds)) {
        dataset <- dataset[, object@sc3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    markers <- object@sc3$biology[[as.character(k)]]$markers
    
    markers <- markers[markers$AUC >= auroc & markers$p.value < p.val, ]
    
    mark.res.plot <- mark_gene_heatmap_param(markers)
    
    row.ann <- data.frame(Cluster = factor(mark.res.plot$sc3_clusters, levels = unique(mark.res.plot$sc3_clusters)))
    rownames(row.ann) <- rownames(mark.res.plot)
    
    do.call(pheatmap::pheatmap, c(list(dataset[rownames(mark.res.plot), , drop = FALSE], 
        show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = res$hc, cutree_cols = k, 
        annotation_row = row.ann, annotation_names_row = FALSE, gaps_row = which(diff(mark.res.plot$sc3_clusters) != 
            0), cellheight = 10), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname sc3_plot_markers.SCESet
#' @aliases sc3_plot_markers
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_plot_markers", signature(object = "SCESet"), function(object, k, auroc = 0.85, 
    p.val = 0.01, show_pdata = NULL) {
    sc3_plot_markers.SCESet(object, k, auroc, p.val, show_pdata)
})

#' Plot stability of the clusters
#' 
#' Stability index shows how stable each cluster is accross the selected 
#' range of ks. The stability index varies between 0 and 1, where 1 means that 
#' the same cluster appears in every solution for different k.
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' 
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw labs ylim
#' 
#' @export
sc3_plot_cluster_stability.SCESet <- function(object, k) {
    
    if (!as.character(k) %in% names(object@sc3$consensus)) {
        stop(paste0("Please choose k from: ", paste(names(object@sc3$consensus), 
            collapse = " ")))
    }
    
    # calculate stability of the clusters check if there are more than 1 k value in
    # ks range
    stability <- NULL
    stability <- StabilityIndex(object@sc3$consensus, k)
    
    d <- data.frame(Cluster = factor(1:length(stability)), Stability = stability)
    ggplot(d, aes(x = d$Cluster, y = d$Stability)) + geom_bar(stat = "identity") + 
        ylim(0, 1) + labs(x = "Cluster", y = "Stability Index") + theme_bw()
}

#' @rdname sc3_plot_cluster_stability.SCESet
#' @aliases sc3_plot_cluster_stability
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_plot_cluster_stability", signature(object = "SCESet"), function(object, 
    k) {
    sc3_plot_cluster_stability.SCESet(object, k)
})
