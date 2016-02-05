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

#' @importFrom robustbase covMcd
#' @importFrom rrcov PcaHubert
#' @importFrom stats qchisq
outl_cells_main <- function(d, chisq.quantile) {
    outl.res <- list()
    for(i in unique(colnames(d))) {
        # reduce p dimensions by using robust PCA
        t <- tryCatch({
            PcaHubert(d[ , colnames(d) == i])
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
                out <- rep(0, dim(d[ , colnames(d) == i])[2])
                names(out) <- rep(i, dim(d[ , colnames(d) == i])[2])
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
                    out <- rep(0, dim(d[ , colnames(d) == i])[2])
                    names(out) <- rep(i, dim(d[ , colnames(d) == i])[2])
                    outl.res[[i]] <- out
                }
            }
        } else {
            out <- rep(0, dim(d[ , colnames(d) == i])[2])
            names(out) <- rep(i, dim(d[ , colnames(d) == i])[2])
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

#' @importFrom stats p.adjust
get_marker_genes <- function(dataset, labels, auroc.threshold) {
    geneAUCs <- apply(dataset, 1, getAUC, labels = labels)
    geneAUCsdf <- data.frame(matrix(unlist(geneAUCs), nrow=length(geneAUCs)/3,
                                    byrow=TRUE))
    rownames(geneAUCsdf) <- rownames(dataset)
    colnames(geneAUCsdf) <- c("AUC","clusts", "p.value")
    geneAUCsdf$AUC <- as.numeric(as.character(geneAUCsdf$AUC))
    geneAUCsdf$clusts <- as.numeric(as.character(geneAUCsdf$clusts))
    geneAUCsdf$p.value <- as.numeric(as.character(geneAUCsdf$p.value))

    geneAUCsdf$p.value <- p.adjust(geneAUCsdf$p.value)
    geneAUCsdf <-
        geneAUCsdf[geneAUCsdf$p.value < 0.01 & !is.na(geneAUCsdf$p.value), ]

    geneAUCsdf <- geneAUCsdf[geneAUCsdf$AUC > auroc.threshold, ]

    d <- NULL
    for(i in sort(unique(geneAUCsdf$clusts))) {
        tmp <- geneAUCsdf[geneAUCsdf$clusts == i, ]
        tmp <- tmp[order(tmp$AUC, decreasing = TRUE),]
        d <- rbind(d, tmp)
    }

    return(d)
}

#' @importFrom stats kruskal.test p.adjust
kruskal_statistics <- function(dataset, labels) {
    t <- apply(dataset, 1, kruskal.test, g = factor(labels))
    ps <- unlist(lapply(t, "[[", "p.value"))
    ps <- p.adjust(ps)
    ps <- ps[!is.na(ps)]
    ps <- ps[ps < 0.01]
    return(ps[order(ps)])
}

