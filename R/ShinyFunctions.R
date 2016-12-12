#' Reindex cluster labels in ascending order
#' 
#' Given an \code{\link[stats]{hclust}} object and the number of clusters \code{k}
#' this function reindex the clusters inferred by \code{cutree(hc, k)[hc$order]}, so that
#' they appear in ascending order. This is particularly useful when plotting
#' heatmaps in which the clusters should be numbered from left to right.
#' 
#' @param hc an object of class hclust
#' @param k number of cluster to be inferred from hc
#'
#' @importFrom stats cutree
#' 
#' @examples
#' hc <- hclust(dist(USArrests), 'ave')
#' cutree(hc, 10)[hc$order]
#' reindex_clusters(hc, 10)[hc$order]
#' 
#' @export
reindex_clusters <- function(hc, k) {
    clusts <- stats::cutree(hc, k)
    labels <- names(clusts)
    names(clusts) <- 1:length(clusts)
    ordering <- clusts[hc$order]
    new.index <- NULL
    j <- 1
    for (i in unique(ordering)) {
        tmp <- rep(j, length(ordering[ordering == i]))
        names(tmp) <- names(ordering[ordering == i])
        new.index <- c(new.index, tmp)
        j <- j + 1
    }
    clusts <- new.index
    clusts <- clusts[order(as.numeric(names(clusts)))]
    names(clusts) <- labels
    return(clusts)
}

#' Find cell outliers in each cluster.
#'
#' Outlier cells in each cluster are detected using robust distances, calculated 
#' using the minimum covariance determinant (MCD), namely using 
#' \code{\link[robustbase]{covMcd}}. The outlier score shows how 
#' different a cell is from all other cells in the cluster and it is defined as 
#' the differences between the square root of the robust distance and the square 
#' root of the 99.99% quantile of the Chi-squared distribution. 
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix
#' @return a numeric vector containing the cell labels and 
#' correspoding outlier scores ordered by the labels
#' @examples
#' d <- get_outl_cells(treutlein[1:10,], colnames(treutlein))
#' head(d)
#' 
#' @importFrom robustbase covMcd
#' @importFrom rrcov PcaHubert
#' @importFrom stats qchisq
#' 
#' @export
get_outl_cells <- function(dataset, labels) {
    chisq.quantile <- 0.9999
    out <- rep(0, length(labels))
    for (i in sort(as.numeric(unique(labels)))) {
        n.cells <- length(labels[labels == i])
        # reduce p dimensions by using robust PCA
        t <- tryCatch({
            PcaHubert(dataset[, labels == i])
        }, warning = function(cond) {
            message(cond)
        }, error = function(cond) {
            message(paste0("No outliers detected in cluster ", i, ". Distribution of gene expression in cells is too skewed towards 0."))
            return(NULL)
        })
        if (class(t) != "NULL") {
            # degrees of freedom used in mcd and chisquare distribution
            if (dim(t@loadings)[1] <= 6) {
                message(paste0("No outliers detected in cluster ", i, ". Small number of cells in the cluster."))
            } else {
                df <- ifelse(dim(t@loadings)[2] > 3, 3, dim(t@loadings)[2])
                
                mcd <- NULL
                if (df != 1) {
                  mcd <- tryCatch({
                    covMcd(t@loadings[, 1:df])
                  }, warning = function(cond) {
                    message(cond)
                  }, error = function(cond) {
                    message("No outliers detected in the cluster. Error in MCD.")
                    return(NULL)
                  })
                }
                
                if (class(mcd) != "NULL") {
                  # sqrt(mcd$mah) - sqrt of robust distance sqrt(qchisq(.95, df = length(mcd$best))) - sqrt of
                  # 97.5% quantile of a chi-squared distribution with p degrees of freedom
                  outliers <- sqrt(mcd$mah) - sqrt(qchisq(chisq.quantile, df = df))
                  outliers[which(outliers < 0)] <- 0
                  out[labels == i] <- outliers
                }
            }
        }
    }
    
    return(out)
}

#' Calculate the area under the ROC curve for a given gene.
#' 
#' For a given gene a binary classifier is constructed 
#' based on the mean cluster expression values (these are calculated using the
#' cell labels). The classifier prediction 
#' is then calculated using the gene expression ranks. The area under the 
#' receiver operating characteristic (ROC) curve is used to quantify the accuracy 
#' of the prediction. A \code{p-value} is assigned to each gene by using the Wilcoxon 
#' signed rank test. 
#' 
#' @param gene expression data of a given gene
#' @param labels cell labels correspodning to the expression values of the gene
#' 
#' @importFrom ROCR prediction performance
#' @importFrom stats aggregate wilcox.test
get_auroc <- function(gene, labels) {
    score <- rank(gene)
    # Get average score for each cluster
    ms <- aggregate(score ~ labels, FUN = mean)
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return NAs if there is a tie for cluster with highest average score (by definition this is
    # not cluster specific)
    if (length(posgroup) > 1) {
        return(c(NA, NA, NA))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest average score vs
    # everything else
    truth <- as.numeric(labels == posgroup)
    # Make predictions & get auc using RCOR package.
    pred <- prediction(score, truth)
    val <- unlist(performance(pred, "auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1], score[truth == 0])$p.value)
    return(c(val, posgroup, pval))
}

#' Calculate marker genes
#'
#' Find marker genes in the dataset. The \code{\link{get_auroc}} is used to calculate
#' marker values for each gene.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding clusters
#' @return data.frame containing the marker genes, corresponding cluster indexes
#' and adjusted \code{p-value}s
#' @importFrom stats p.adjust
#' @examples
#' d <- get_marker_genes(treutlein[1:10,], colnames(treutlein))
#' d
#' 
#' @export
get_marker_genes <- function(dataset, labels) {
    res <- apply(dataset, 1, get_auroc, labels = labels)
    res <- data.frame(matrix(unlist(res), ncol = 3, byrow = T))
    colnames(res) <- c("auroc", "clusts", "pvalue")
    res$pvalue <- p.adjust(res$pvalue)
    return(res)
}

#' Get marker genes from an object of \code{SCESet} class
#' 
#' This functions returns all marker gene columns from the \code{phenoData} slot 
#' of the input object corresponding to the number of clusters \code{k}. Additionally,
#' it rearranges genes by the cluster index and order them by the area under the 
#' ROC curve value inside of each cluster.
#'
#' @param object an objec of \code{SCESet} class
#' @param k number of cluster
#' @param p_val p-value threshold
#' @param auroc area under the ROC curve threshold
#'
organise_marker_genes <- function(object, k, p_val, auroc) {
    dat <- object@featureData@data[, c(paste0("sc3_", k, "_markers_clusts"), paste0("sc3_", k, 
        "_markers_auroc"), paste0("sc3_", k, "_markers_padj"))]
    dat <- dat[dat[, paste0("sc3_", k, "_markers_padj")] < p_val & !is.na(dat[, paste0("sc3_", 
        k, "_markers_padj")]), ]
    dat <- dat[dat[, paste0("sc3_", k, "_markers_auroc")] > auroc, ]
    
    d <- NULL
    
    for (i in sort(unique(dat[, paste0("sc3_", k, "_markers_clusts")]))) {
        tmp <- dat[dat[, paste0("sc3_", k, "_markers_clusts")] == i, ]
        tmp <- tmp[order(tmp[, paste0("sc3_", k, "_markers_auroc")], decreasing = TRUE), ]
        d <- rbind(d, tmp)
    }
    
    return(d)
}

#' Reorder and subset gene markers for plotting on a heatmap
#' 
#' Reorders the rows of the input data.frame based on the \code{sc3_k_markers_clusts}
#' column and also keeps only the top 10 genes for each value of \code{sc3_k_markers_clusts}.
#'
#' @param markers a \code{data.frame} object with the following colnames:
#' \code{sc3_k_markers_clusts}, \code{sc3_k_markers_auroc}, \code{sc3_k_markers_padj}.
#' 
markers_for_heatmap <- function(markers) {
    res <- NULL
    for (i in unique(markers[, 1])) {
        tmp <- markers[markers[, 1] == i, ]
        if (nrow(tmp) > 10) {
            res <- rbind(res, tmp[1:10, ])
        } else {
            res <- rbind(res, tmp)
        }
    }
    
    return(res)
}

#' Find differentially expressed genes
#'
#' Differential expression is calculated using the non-parametric Kruskal-Wallis test.
#' A significant \code{p-value} indicates that gene expression in at least one cluster
#' stochastically dominates one other cluster. Note that the calculation of 
#' differential expression after clustering can introduce a bias in the 
#' distribution of \code{p-value}s, and thus we advise to use the \code{p-value}s 
#' for ranking the genes only.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix
#' @return a numeric vector containing the differentially expressed genes and 
#' correspoding p-values
#' @examples
#' d <- get_de_genes(treutlein[1:10, ], colnames(treutlein))
#' head(d)
#' 
#' @importFrom stats kruskal.test p.adjust
#' 
#' @export
get_de_genes <- function(dataset, labels) {
    tmp <- apply(dataset, 1, kruskal.test, g = factor(labels))
    ps <- unlist(lapply(tmp, "[[", "p.value"))
    ps <- p.adjust(ps)
    return(ps)
}

#' Get differentiall expressed genes from an object of \code{SCESet} class
#' 
#' This functions returns all marker gene columns from the \code{phenoData} slot 
#' of the input object corresponding to the number of clusters \code{k}. Additionally,
#' it rearranges genes by the cluster index and order them by the area under the 
#' ROC curve value inside of each cluster.
#'
#' @param object an objec of \code{SCESet} class
#' @param k number of cluster
#' @param p_val p-value threshold
#' 
organise_de_genes <- function(object, k, p_val) {
    de_genes <- object@featureData@data[, paste0("sc3_", k, "_de_padj")]
    names(de_genes) <- rownames(object@featureData@data)
    de_genes <- de_genes[!is.na(de_genes)]
    de_genes <- de_genes[de_genes < p_val]
    de_genes <- de_genes[order(de_genes)]
    return(de_genes)
}

#' Calculate the stability index of the obtained clusters when changing k
#'
#' Stability index shows how stable each cluster is accross the selected range of k. 
#' The stability index varies between 0 and 1, where 
#' 1 means that the same cluster appears in every solution for different k.
#' 
#' Formula (imagine a given cluster with is split into N clusters when k is changed, and
#' in each of the new clusters there are given_cells of the given cluster and also some extra_cells from other clusters):
#' SI = sum_over_ks(sum_over_clusters_N(given_cells/(given_cells + extra_cells)))/N(corrects for stability of each cluster)/N(corrects for the number of clusters)/length(ks)
#'
#' @param consensus consensus item of the sc3 slot of an object of 'SCESet' class
#' @param k number of clusters k
#' @return a numeric vector containing a stability index of each cluster
calculate_stability <- function(consensus, k) {
    hc <- consensus[[as.character(k)]]$hc
    labs <- reindex_clusters(hc, k)
    
    kMax <- max(as.numeric(names(consensus)))
    kMin <- min(as.numeric(names(consensus)))
    kRange <- kMax - kMin
    
    stability <- rep(0, k)
    
    for (i in 1:k) {
        inds <- names(labs[labs == i])
        # sum over k range
        for (k2 in kMin:kMax) {
            if (k2 != k) {
                hc2 <- consensus[[as.character(k2)]]$hc
                labs2 <- cutree(hc2, k2)
                clusts <- as.numeric(names(table(labs2[names(labs2) %in% inds])))
                N <- length(clusts)
                # sum over new clusters, taking into account new cells from other clusters
                for (j in clusts) {
                  inds2 <- names(labs2[labs2 == j])
                  s <- length(inds[inds %in% inds2])/length(inds2)/N/N/kRange
                  stability[i] <- stability[i] + s
                }
            }
        }
    }
    return(stability)
}
