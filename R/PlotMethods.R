sc3_plot_consensus <- function(object, k) {
    
    res <- prepare_output(object, k)
    
    consensus <- object@consensus$sc3_consensus[[as.character(k)]]$consensus
    
    pheatmap::pheatmap(
        consensus,
        cluster_rows = res$hc,
        cluster_cols = res$hc,
        cutree_rows = k,
        cutree_cols = k,
        show_rownames = FALSE,
        show_colnames = FALSE
    )
}

sc3_plot_expression <- function(object, k) {
    
    res <- prepare_output(object, k)
    
    dataset <- object@consensus$sc3_processed_dataset
    
    if(!is.null(object@consensus$svm_train_inds)) {
        dataset <- dataset[ , object@consensus$svm_train_inds]
    }
    
    pheatmap::pheatmap(
        dataset,
        cluster_cols = res$hc,
        kmeans_k = 100,
        cutree_cols = k,
        show_rownames = FALSE,
        show_colnames = FALSE
    )
}

sc3_plot_tsne <- function(object, k, perplexity = floor(ncol(object@consensus$sc3_processed_dataset) / 5), seed = 1234567) {
    res <- prepare_output(object, k)
    
    dataset <- object@consensus$sc3_processed_dataset
    
    if(!is.null(object@consensus$svm_train_inds)) {
        dataset <- dataset[ , object@consensus$svm_train_inds]
    }
    
    dataset <- dataset[ , res$hc$order]
    colnames(dataset) <- res$new.labels
    
    set.seed(seed)
    tsne_out <- Rtsne::Rtsne(
        t(dataset),
        perplexity = perplexity
    )
    df_to_plot <- as.data.frame(tsne_out$Y)
    df_to_plot$Cluster <- factor(
        colnames(dataset),
        levels = unique(colnames(dataset))
    )
    comps <- colnames(df_to_plot)[1:2]
    ggplot(df_to_plot, aes_string(x = comps[1],
                                  y = comps[2],
                                  color = "Cluster")) +
        geom_point() +
        xlab("Dimension 1") +
        ylab("Dimension 2") +
        theme_bw()
}

sc3_plot_de_genes <- function(object, k, p.val = 0.01) {
    
    res <- prepare_output(object, k)
    
    dataset <- object@consensus$sc3_processed_dataset
    
    if(!is.null(object@consensus$svm_train_inds)) {
        dataset <- dataset[ , object@consensus$svm_train_inds]
    }
    
    de.genes <- object@consensus$sc3_biology[[as.character(k)]]$de.genes
    
    de.genes <- de.genes[de.genes$p.value < p.val, , drop = FALSE]
    d <- head(de.genes, 50)
    row.ann <- data.frame("p.value" = -log10(d))
    rownames(row.ann) <- rownames(d)
    pheatmap::pheatmap(
        dataset[rownames(d), , drop = FALSE],
        show_colnames = FALSE,
        cluster_rows = FALSE,
        cluster_cols = res$hc,
        cutree_cols = k,
        annotation_row = row.ann,
        annotation_names_row = FALSE,
        cellheight = 10
    )
}

sc3_plot_markers <- function(object, k, auroc = 0.85, p.val = 0.01) {
    
    res <- prepare_output(object, k)
    
    dataset <- object@consensus$sc3_processed_dataset
    
    if(!is.null(object@consensus$svm_train_inds)) {
        dataset <- dataset[ , object@consensus$svm_train_inds]
    }
    
    markers <- object@consensus$sc3_biology[[as.character(k)]]$markers
    
    markers <- markers[markers$AUC >= auroc & markers$p.value < p.val, ]

    mark.res.plot <- mark_gene_heatmap_param(markers)
    
    row.ann <- data.frame(
        Cluster = factor(
            mark.res.plot$clusts,
            levels = unique(mark.res.plot$clusts)
        )
    )
    rownames(row.ann) <- rownames(mark.res.plot)
    
    pheatmap::pheatmap(
        dataset[rownames(mark.res.plot), , drop = FALSE],
        show_colnames = FALSE,
        cluster_rows = FALSE,
        cluster_cols = res$hc,
        cutree_cols = k,
        annotation_row = row.ann,
        annotation_names_row = FALSE,
        gaps_row = which(diff(mark.res.plot$clusts) != 0),
        cellheight = 10
    )
}

sc3_plot_cell_outliers <- function(object, k) {
    res <- prepare_output(object, k)
    
    outl <- object@consensus$sc3_biology[[as.character(k)]]$cell.outl
    outl <- outl[res$hc$order, ]
    
    outl$cell.ind <- 1:nrow(outl)
    
    cols <- iwanthue(length(unique(outl$clusts)))
    
    outl$clusts <- factor(
        outl$clusts,
        levels =
            unique(
                as.character(
                    outl$clusts
                )
            )
    )
    
    comps <- colnames(outl)
    
    ggplot(outl, aes_string(x = comps[3],
                  y = comps[2],
                  fill = comps[1], 
                  color = comps[1])) +
        geom_bar(stat = "identity") +
        geom_point() +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        guides(color = FALSE, fill = FALSE) +
        labs(x = "Cells", y = "Outlier score") +
        theme_bw()
}

sc3_plot_cluster_stability <- function(object, k) {
    
    if(!as.character(k) %in% names(object@consensus$sc3_consensus)) {
        stop(paste0("Please choose k from: ", paste(names(object@consensus$sc3_consensus), collapse = " "))) 
    }
    
    res <- prepare_output(object, k)
    
    # calculate stability of the clusters
    # check if there are more than 1 k value in ks range
    stability <- NULL
    stability <- StabilityIndex(
        object,
        k
    )
    
    d <- data.frame(
        Cluster = factor(1:length(stability)),
        Stability = stability)
    ggplot(d, aes(x = d$Cluster, y = d$Stability)) +
        geom_bar(stat = "identity") +
        ylim(0, 1) +
        labs(x = "Cluster", y = "Stability Index") +
        theme_bw()
}

