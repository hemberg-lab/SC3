reindex_clusters <- function(ordering) {
    new.index <- NULL
    j <- 1
    for(i in unique(ordering)) {
        tmp <- rep(j, length(ordering[ordering == i]))
        names(tmp) <- names(ordering[ordering == i])
        new.index <- c(new.index, tmp)
        j <- j + 1
    }
    return(new.index)
}

de_gene_heatmap_param <- function(res) {
    row.ann <- data.frame("minus.log10.p.value" = -log10(res))
    rownames(row.ann) <- names(res)

    return(list(row.ann = row.ann))
}

mark_gene_heatmap_param <- function(mark.res, labs) {
    mark.res.plot <- NULL
    for(i in labs) {
        tmp <- mark.res[mark.res[,2] == i, ]
        if(dim(tmp)[1] > 10) {
            mark.res.plot <- rbind(mark.res.plot, tmp[1:10, ])
        } else {
            mark.res.plot <- rbind(mark.res.plot, tmp)
        }
    }

    row.ann <- data.frame(Cluster =
                              factor(mark.res.plot$clusts,levels =
                                         unique(mark.res.plot$clusts)))
    rownames(row.ann) <- rownames(mark.res.plot)

    row.gaps <- as.numeric(mark.res.plot$clusts)
    row.gaps <- which(diff(row.gaps) != 0)

    return(list(mark.res.plot = mark.res.plot, row.ann = row.ann,
                row.gaps = row.gaps))
}

#' Find cell outliers
#'
#' If the cell labels are available this functions allows a user to calculate
#' cell outlier scores manually.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix
#' @param chisq.quantile a threshold of the chi-squared distribution used for 
#' cell outliers detection, default is 0.9999
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
get_outl_cells <- function(dataset, labels, chisq.quantile = 0.9999) {
    outl.res <- list()
    for(i in unique(labels)) {
        n.cells <- length(labels[labels == i])
        # reduce p dimensions by using robust PCA
        t <- tryCatch({
            PcaHubert(dataset[ , labels == i])
        }, warning = function(cond) {
            message(cond)
        }, error = function(cond) {
            message(paste0("No outliers detected in cluster ", i,
                           ". Distribution of gene expression in cells is too skewed towards 0."))
            return(NULL)
        })
        if(class(t) != "NULL") {
            # degrees of freedom used in mcd and chisquare distribution
            if(dim(t@loadings)[1] <= 6) {
                message(paste0("No outliers detected in cluster ", i,
                               ". Small number of cells in the cluster."))
                out <- rep(0, n.cells)
                names(out) <- rep(i, n.cells)
                outl.res[[i]] <- out
            } else {
                df <- ifelse(dim(t@loadings)[2] > 3, 3, dim(t@loadings)[2])

                mcd <- NULL
                if(df != 1) {
                    mcd <- tryCatch({
                        covMcd(t@loadings[ , 1:df])
                    }, warning = function(cond) {
                        message(cond)
                    }, error = function(cond) {
                        message("No outliers detected in the cluster. Error in MCD.")
                        return(NULL)
                    })
                }

                if(class(mcd) != "NULL") {
                    # sqrt(mcd$mah) - sqrt of robust distance
                    # sqrt(qchisq(.95, df = length(mcd$best))) -
                    # sqrt of 97.5% quantile of a
                    # chi-squared distribution with p degrees of freedom
                    outliers <-
                        sqrt(mcd$mah) - sqrt(qchisq(chisq.quantile, df = df))
                    outliers[which(outliers < 0)] <- 0
                    outl.res[[i]] <- outliers
                } else {
                    out <- rep(0, n.cells)
                    names(out) <- rep(i, n.cells)
                    outl.res[[i]] <- out
                }
            }
        } else {
            out <- rep(0, n.cells)
            names(out) <- rep(i, n.cells)
            outl.res[[i]] <- out
        }
    }

    nams <- NULL
    vals <- NULL
    for(i in 1:length(outl.res)) {
        vals <- c(vals, outl.res[[i]])
        nams <- c(nams, names(outl.res[[i]]))
    }
    names(vals) <- nams

    return(vals)
}

#' @importFrom RSelenium remoteDriver
open_webgestalt_go <- function(genes) {
    remDr <- remoteDriver(remoteServerAddr = "localhost"
                          , port = 4444
                          , browserName = "firefox"
    )
    remDr$open()
    remDr$navigate("http://bioinfo.vanderbilt.edu/webgestalt/login.php")
    webElem <- remDr$findElement(using = 'id', value = "email")
    webElem$sendKeysToElement(list("vk6@sanger.ac.uk"))
    webElem <- remDr$findElement(using = 'name', "remember")
    webElem$clickElement()
    webElem <- remDr$findElement(using = 'name', "submit")
    webElem$clickElement()
    webElem <- remDr$findElement(using = 'name', "pastefile")
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
    # Return negatives if there is a tie for cluster with highest average score
    # (by definition this is not cluster specific)
    if(length(posgroup) > 1) {
        return (c(-1,-1,1))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest
    # average score vs everything else
    truth <- as.numeric(labels == posgroup)
    #Make predictions & get auc using RCOR package.
    pred <- prediction(score,truth)
    val <- unlist(performance(pred,"auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1],
                                         score[truth == 0])$p.value)
    return(c(val,posgroup,pval))
}

#' Find marker genes
#'
#' If the cell labels are available this functions allows a user to calculate
#' marker genes manually.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix.
#' Labels must be integers or character integers, e.g. 1, 2, 3 or "1", "2", "3" ect.
#' @param auroc.threshold area under the ROC curve threshold, by default it is
#' 0.85. Values close to 0.5 will include very weak marker genes, values close
#' to 1 will only include very strong marker genes.
#' @param p.val p-value threshold, by default it is 0.01
#' @return data.frame containing the marker genes
#' @importFrom stats p.adjust
#' @examples
#' d <- get_marker_genes(treutlein, colnames(treutlein))
#' head(d)
#' 
#' @export
get_marker_genes <- function(dataset, labels, auroc.threshold = 0.85, p.val = 0.01) {
    geneAUCs <- apply(dataset, 1, getAUC, labels = labels)
    geneAUCsdf <- data.frame(matrix(unlist(geneAUCs), nrow=length(geneAUCs)/3,
                                    byrow=TRUE))
    rownames(geneAUCsdf) <- rownames(dataset)
    colnames(geneAUCsdf) <- c("AUC","clusts", "p.value")
    # remove genes with ties
    geneAUCsdf <- geneAUCsdf[geneAUCsdf$clusts != -1, ]
    geneAUCsdf$AUC <- as.numeric(as.character(geneAUCsdf$AUC))
    geneAUCsdf$clusts <- as.numeric(as.character(geneAUCsdf$clusts))
    geneAUCsdf$p.value <- as.numeric(as.character(geneAUCsdf$p.value))

    geneAUCsdf$p.value <- p.adjust(geneAUCsdf$p.value)
    geneAUCsdf <-
        geneAUCsdf[geneAUCsdf$p.value < p.val & !is.na(geneAUCsdf$p.value), ]

    geneAUCsdf <- geneAUCsdf[geneAUCsdf$AUC > auroc.threshold, ]

    d <- NULL
    for(i in sort(unique(geneAUCsdf$clusts))) {
        tmp <- geneAUCsdf[geneAUCsdf$clusts == i, ]
        tmp <- tmp[order(tmp$AUC, decreasing = TRUE),]
        d <- rbind(d, tmp)
    }

    return(d)
}

#' Find differentially expressed genes
#'
#' If the cell labels are available this functions allows a user to calculate
#' differentially expressed genes manually.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding to the columns of the expression matrix
#' @param p.val p-value threshold, by default it is 0.01
#' @return a numeric vector containing the differentially expressed genes and 
#' correspoding p-values
#' @examples
#' d <- get_de_genes(treutlein, colnames(treutlein))
#' head(d)
#' 
#' @importFrom stats kruskal.test p.adjust
#' 
#' @export
get_de_genes <- function(dataset, labels, p.val = 0.01) {
    t <- apply(dataset, 1, kruskal.test, g = factor(labels))
    ps <- unlist(lapply(t, "[[", "p.value"))
    ps <- p.adjust(ps)
    ps <- ps[!is.na(ps)]
    ps <- ps[ps < p.val]
    return(ps[order(ps)])
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
#' @param stab.res internal matrix of precomputed clustering results
#' @param k current value of the number of clusters k
#' @return a numeric vector containing a stability index of each cluster
StabilityIndex <- function(stab.res, k) {
    hc <- stab.res[as.numeric(stab.res[ , 1]) == k, 2][[1]][[3]]
    labs <- cutree(hc, k)
    labs <- labs[hc$order]
    labs <- reindex_clusters(labs)
    
    kMax <- max(unique(as.numeric(stab.res[ , 1])))
    kMin <- min(unique(as.numeric(stab.res[ , 1])))
    kRange <- kMax - kMin
    
    stability <- rep(0, k)
    
    for (i in 1:k) {
        inds <- names(labs[labs == i])
        # sum over k range
        for (k2 in kMin:kMax) {
            if (k2 != k) {
                hc2 <- stab.res[as.numeric(stab.res[ , 1]) == k2, 2][[1]][[3]]
                labs2 <- cutree(hc2, k2)
                clusts <- as.numeric(names(table(labs2[names(labs2) %in% inds])))
                N <- length(clusts)
                # sum over new clusters, taking into account new cells from
                # other clusters
                for(j in clusts) {
                    inds2 <- names(labs2[labs2 == j])
                    s <- length(inds[inds %in% inds2]) / length(inds2) / N / N / kRange
                    stability[i] <- stability[i] + s
                }
            }
        }
    }
    return(stability)
}
