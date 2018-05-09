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
#' @name sc3_plot_consensus
#' @aliases sc3_plot_consensus, sc3_plot_consensus,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
sc3_plot_consensus.SingleCellExperiment <- function(object, k, show_pdata) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
    consensus <- metadata(object)$sc3$consensus[[as.character(k)]]$consensus
    
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
    do.call(pheatmap::pheatmap, c(list(consensus, cluster_rows = hc, cluster_cols = hc, cutree_rows = k, 
        cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname sc3_plot_consensus
#' @aliases sc3_plot_consensus
setMethod("sc3_plot_consensus", signature(object = "SingleCellExperiment"), sc3_plot_consensus.SingleCellExperiment)

#' Plot silhouette indexes of the cells
#' 
#' A silhouette is a quantitative measure of the diagonality of the consensus 
#' matrix. An average silhouette width (shown at the bottom left of the silhouette 
#' plot) varies from 0 to 1, where 1 represents a perfectly block-diagonal 
#' consensus matrix and 0 represents a situation where there is no 
#' block-diagonal structure. The best clustering is achieved when the average 
#' silhouette width is close to 1.
#' 
#' @name sc3_plot_silhouette
#' @aliases sc3_plot_silhouette, sc3_plot_silhouette,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
sc3_plot_silhouette.SingleCellExperiment <- function(object, k) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    silh <- metadata(object)$sc3$consensus[[as.character(k)]]$silhouette
    plot(silh, col = "black")
}

#' @rdname sc3_plot_silhouette
#' @aliases sc3_plot_silhouette
setMethod("sc3_plot_silhouette", signature(object = "SingleCellExperiment"), sc3_plot_silhouette.SingleCellExperiment)

#' Plot expression matrix used for SC3 clustering as a heatmap
#' 
#' The expression panel represents the original input expression matrix 
#' (cells in columns and genes in rows) after the gene filter. 
#' Genes are clustered by kmeans with k = 100 (dendrogram on the left) and 
#' the heatmap represents the expression levels of the gene cluster centers 
#' after log2-scaling.
#' 
#' @name sc3_plot_expression
#' @aliases sc3_plot_expression, sc3_plot_expression,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
sc3_plot_expression.SingleCellExperiment <- function(object, k, show_pdata) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$sc3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
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
    
    if(nrow(dataset) > 100) {
        do.call(pheatmap::pheatmap, c(list(dataset, cluster_cols = hc, kmeans_k = 100, cutree_cols = k, 
                                           show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
    } else {
        do.call(pheatmap::pheatmap, c(list(dataset, cluster_cols = hc, cutree_cols = k, 
                                           show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
    }
}

#' @rdname sc3_plot_expression
#' @aliases sc3_plot_expression
setMethod("sc3_plot_expression", signature(object = "SingleCellExperiment"), sc3_plot_expression.SingleCellExperiment)

#' Plot expression of DE genes of the clusters identified by \code{SC3} as a heatmap
#' 
#' \code{SC3} plots gene expression profiles of the 50 genes with the lowest p-values. 
#' 
#' @name sc3_plot_de_genes
#' @aliases sc3_plot_de_genes, sc3_plot_de_genes,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
sc3_plot_de_genes.SingleCellExperiment <- function(object, k, p.val, show_pdata) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$sc3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
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
    # remove Inf when the p-value is actually 0 (less than the accuracy limit)
    de_genes[de_genes < 1e-17] <- 1e-17
    row_ann <- data.frame(log10_padj = -log10(de_genes))
    rownames(row_ann) <- names(de_genes)
    
    do.call(pheatmap::pheatmap, c(list(dataset[names(de_genes), , drop = FALSE], show_colnames = FALSE, 
        cluster_rows = FALSE, cluster_cols = hc, cutree_cols = k, annotation_row = row_ann, cellheight = 10), 
        list(annotation_col = ann)[add_ann_col]))
}

#' @rdname sc3_plot_de_genes
#' @aliases sc3_plot_de_genes
setMethod("sc3_plot_de_genes", signature(object = "SingleCellExperiment"), sc3_plot_de_genes.SingleCellExperiment)

#' Plot expression of marker genes identified by \code{SC3} as a heatmap.
#' 
#' By default the genes with the area under the ROC curve (AUROC) > 0.85 
#' and with the p-value < 0.01 are selected and the top 10 marker 
#' genes of each cluster are visualized in this heatmap.
#' 
#' @name sc3_plot_markers
#' @aliases sc3_plot_markers, sc3_plot_markers,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param auroc area under the ROC curve
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
sc3_plot_markers.SingleCellExperiment <- function(object, k, auroc, p.val, show_pdata) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$sc3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
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
    
    # get all marker genes
    markers <- organise_marker_genes(object, k, p.val, auroc)
    
    if(!is.null(markers)) {
        # get top 10 marker genes of each cluster
        markers <- markers_for_heatmap(markers)
        
        row.ann <- data.frame(Cluster = factor(markers[, 1], levels = unique(markers[, 1])))
        rownames(row.ann) <- markers$feature_symbol
        
        do.call(pheatmap::pheatmap, c(list(dataset[markers$feature_symbol, , drop = FALSE], show_colnames = FALSE, 
            cluster_rows = FALSE, cluster_cols = hc, cutree_cols = k, annotation_row = row.ann, annotation_names_row = FALSE, 
            gaps_row = which(diff(markers[, 1]) != 0), cellheight = 10), list(annotation_col = ann)[add_ann_col]))
    } else {
        message("No markers have been found, try to lower significance thresholds!")
    }
}

#' @rdname sc3_plot_markers
#' @aliases sc3_plot_markers
setMethod("sc3_plot_markers", signature(object = "SingleCellExperiment"), sc3_plot_markers.SingleCellExperiment)

#' Plot stability of the clusters
#' 
#' Stability index shows how stable each cluster is accross the selected 
#' range of ks. The stability index varies between 0 and 1, where 1 means that 
#' the same cluster appears in every solution for different k.
#' 
#' @name sc3_plot_cluster_stability
#' @aliases sc3_plot_cluster_stability, sc3_plot_cluster_stability,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' 
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw labs ylim
sc3_plot_cluster_stability.SingleCellExperiment <- function(object, k) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    # calculate stability of the clusters check if there are more than 1 k value in ks range
    stability <- NULL
    stability <- calculate_stability(metadata(object)$sc3$consensus, k)
    
    d <- data.frame(Cluster = factor(1:length(stability)), Stability = stability)
    ggplot(d, aes(x = d$Cluster, y = d$Stability)) + geom_bar(stat = "identity") + ylim(0, 1) + 
        labs(x = "Cluster", y = "Stability Index") + theme_bw()
}

#' @rdname sc3_plot_cluster_stability
#' @aliases sc3_plot_cluster_stability
setMethod("sc3_plot_cluster_stability", signature(object = "SingleCellExperiment"), sc3_plot_cluster_stability.SingleCellExperiment)
