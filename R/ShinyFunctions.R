
#' @importFrom stats cutree
prepare_output <- function(object, k) {
    consensus <- object@sc3$consensus
    ks <- as.numeric(names(consensus))
    # get all results for k
    res <- consensus[[as.character(k)]]
    hc <- res$hc
    silh <- res$silhouette
    
    return(list(hc = hc, silh = silh, cell.order = hc$order))
}

#' @importFrom stats cutree
get_clusts <- function(hc, k) {
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

mark_gene_heatmap_param <- function(markers) {
    mark.res.plot <- NULL
    for (i in unique(markers$sc3_clusters)) {
        tmp <- markers[markers$sc3_clusters == i, ]
        if (nrow(tmp) > 10) {
            mark.res.plot <- rbind(mark.res.plot, tmp[1:10, ])
        } else {
            mark.res.plot <- rbind(mark.res.plot, tmp)
        }
    }
    
    return(mark.res.plot)
}

#' Find cell outliers
#'
#' If the cell labels are available this functions allows a user to calculate
#' cell outlier scores manually.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix
#' @return a numeric vector containing the cell labels and 
#' correspoding outlier scores ordered by the labels
#' @examples
#' d <- get_outl_cells(treutlein, colnames(treutlein))
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
                  # sqrt(mcd$mah) - sqrt of robust distance sqrt(qchisq(.95, df =
                  # length(mcd$best))) - sqrt of 97.5% quantile of a chi-squared distribution with
                  # p degrees of freedom
                  outliers <- sqrt(mcd$mah) - sqrt(qchisq(chisq.quantile, df = df))
                  outliers[which(outliers < 0)] <- 0
                  out[labels == i] <- outliers
                }
            }
        }
    }
    
    return(out)
}

#' @importFrom RSelenium remoteDriver
open_gprofiler <- function(genes) {
    remDr <- remoteDriver(remoteServerAddr = "localhost", port = 4444, browserName = "firefox")
    remDr$open()
    remDr$navigate("http://biit.cs.ut.ee/gprofiler")
    webElem <- remDr$findElement(using = "id", "query")
    webElem$sendKeysToElement(as.list(paste0(genes, "\n")))
}


#' @importFrom ROCR prediction performance
#' @importFrom stats aggregate wilcox.test
getAUC <- function(gene, labels) {
    score <- rank(gene)
    # Get average score for each cluster
    ms <- aggregate(score ~ labels, FUN = mean)
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return negatives if there is a tie for cluster with highest average score (by
    # definition this is not cluster specific)
    if (length(posgroup) > 1) {
        return(c(-1, -1, 1))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest average score
    # vs everything else
    truth <- as.numeric(labels == posgroup)
    # Make predictions & get auc using RCOR package.
    pred <- prediction(score, truth)
    val <- unlist(performance(pred, "auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1], score[truth == 0])$p.value)
    return(c(val, posgroup, pval))
}

#' Find marker genes
#'
#' If the cell labels are available this functions allows a user to calculate
#' marker genes manually.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix.
#' Labels must be integers or character integers, e.g. 1, 2, 3 or '1', '2', '3' ect.
#' @return data.frame containing the marker genes
#' @importFrom stats p.adjust
#' @examples
#' d <- get_marker_genes(treutlein, colnames(treutlein))
#' head(d)
#' 
#' @export
get_marker_genes <- function(dataset, labels) {
    auroc.threshold <- 0.5
    p.val <- 0.1
    geneAUCs <- apply(dataset, 1, getAUC, labels = labels)
    geneAUCsdf <- data.frame(matrix(unlist(geneAUCs), nrow = length(geneAUCs)/3, 
        byrow = TRUE))
    rownames(geneAUCsdf) <- rownames(dataset)
    colnames(geneAUCsdf) <- c("AUC", "sc3_clusters", "p.value")
    # remove genes with ties
    geneAUCsdf <- geneAUCsdf[geneAUCsdf$sc3_clusters != -1, ]
    geneAUCsdf$AUC <- as.numeric(as.character(geneAUCsdf$AUC))
    geneAUCsdf$sc3_clusters <- as.numeric(as.character(geneAUCsdf$sc3_clusters))
    geneAUCsdf$p.value <- as.numeric(as.character(geneAUCsdf$p.value))
    
    geneAUCsdf$p.value <- p.adjust(geneAUCsdf$p.value)
    geneAUCsdf <- geneAUCsdf[geneAUCsdf$p.value < p.val & !is.na(geneAUCsdf$p.value), 
        ]
    
    geneAUCsdf <- geneAUCsdf[geneAUCsdf$AUC > auroc.threshold, ]
    
    d <- NULL
    for (i in sort(unique(geneAUCsdf$sc3_clusters))) {
        tmp <- geneAUCsdf[geneAUCsdf$sc3_clusters == i, ]
        tmp <- tmp[order(tmp$AUC, decreasing = TRUE), ]
        d <- rbind(d, tmp)
    }
    
    return(d)
}

#' Find differentially expressed genes
#'
#' Differential expression is calculated using the non-parametric 
#' Kruskal-Wallis test. A significant p-value indicates that gene 
#' expression in at least one cluster stochastically dominates one other cluster.
#' Note that the calculation of differential 
#' expression after clustering can introduce a bias in the distribution of 
#' p-values, and thus we advise to use the p-values for ranking the genes only.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix
#' @return a numeric vector containing the differentially expressed genes and 
#' correspoding p-values
#' @examples
#' d <- get_de_genes(treutlein, colnames(treutlein))
#' head(d)
#' 
#' @importFrom stats kruskal.test p.adjust
#' 
#' @export
get_de_genes <- function(dataset, labels) {
    t <- apply(dataset, 1, kruskal.test, g = factor(labels))
    ps <- unlist(lapply(t, "[[", "p.value"))
    ps <- p.adjust(ps)
    return(ps)
}

organise_de_genes <- function(object, k, p_val) {
    de_genes <- object@featureData@data[ , paste0("sc3_", k, "_de_padj")]
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
StabilityIndex <- function(consensus, k) {
    hc <- consensus[[as.character(k)]]$hc
    labs <- get_clusts(hc, k)
    
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
