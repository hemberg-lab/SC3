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
    if ( is.null(dataset) ) {
        warning(paste0("The object does not contain ", exprs_values, " expression values. Returning NULL."))
        return(object)
    }
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
sc3_estimate_k.SCESet <- function(
    object
    ) {
    dataset <- object@consensus$sc3_processed_dataset
    if ( is.null(dataset) ) {
        warning(paste0("Please run sc3_process() first!"))
        return(object)
    }
    res <- estkTW(dataset = dataset)
    object@consensus$"sc3_k_prediction" <- res
    return(object)
}

#' @export
setMethod("sc3_estimate_k", signature(object = "SCESet"),
          function(
              object
          ) {
              sc3_estimate_k.SCESet(object)
          })

#' @export
sc3_process.SCESet <- function(
    object,
    exprs_values = "counts",
    gene.filter = TRUE,
    gene.filter.fraction = 0.06,
    gene.reads.rare = 2,
    gene.reads.ubiq = 0,
    log.scale = TRUE,
    d.region.min = 0.04,
    d.region.max = 0.07,
    svm.num.cells = NULL,
    svm.train.inds = NULL,
    n.cores = NULL) {
    dataset <- object@assayData[[exprs_values]]
    if ( is.null(dataset) ) {
        warning(paste0("The object does not contain ", exprs_values, " expression values. Returning NULL."))
        return(object)
    }
    # remove duplicated genes
    message("Removing duplicated genes...")
    dataset <- dataset[!duplicated(rownames(dataset)), ]
    
    # gene filter
    if(gene.filter) {
        dataset <- gene_filter(dataset, gene.filter.fraction, gene.reads.rare, gene.reads.ubiq)
        if(nrow(dataset) == 0) {
            message("All genes were removed after the gene filter! Stopping now...")
            return(object)
        }
    }
    
    # log2 transformation
    if(log.scale) {
        message("log2-scaling...")
        dataset <- log2(1 + dataset)
    }
    
    object@consensus$"sc3_processed_dataset" <- dataset
    
    if(ncol(dataset) > 2000) {
        object@consensus$"sc3_kmeans_nstart" <- 50
        message("Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...")
    } else {
        object@consensus$"sc3_kmeans_nstart" <- 1000
    }
    
    # define number of cells and region of dimensions
    n.dim <- floor(d.region.min * ncol(dataset)) : ceiling(d.region.max * ncol(dataset))
    
    # for large datasets restrict the region of dimensions to 15
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }
    
    object@consensus$"sc3_n_dim" <- n.dim
    
    # prepare for SVM
    if(
        !is.null(svm.num.cells) | 
        !is.null(svm.train.inds) | 
        ncol(dataset) > 5000
    ) {
        # handle all possible errors
        if(!is.null(svm.num.cells)) {
            if(!is.null(svm.train.inds)) {
                return(
                    message(
                        "You have set both svm.num.cells and svm.train.inds parameters for SVM training. Please set only one of them and rerun sc3_process()."
                    )
                )
            }
            if(svm.num.cells >= ncol(dataset) - 1) return(
                message(
                    "Number of cells used for SVM training is larger (or equal) than the total number of cells in your dataset. Please make svm.num.cells parameter smaller and rerun sc3_process()."
                )
            )
            if(svm.num.cells < 10) {
                return(
                    message(
                        "Number of cells used for SVM training is less than 10. Please make sure the number of clusters k is smaller than 10 or increase the number of training cells."
                    )
                )
            }
        }
        if(!is.null(svm.train.inds)) {
            if(length(svm.train.inds) < 10) {
                return(
                    message(
                        "Number of cells used for SVM training is less than 10. Please make sure the number of clusters k is smaller than 10 or increase the number of training cells."
                    )
                )
            }
            if(max(svm.train.inds) > ncol(dataset) - 1) {
                return(
                    message(
                        "Number of cells used for SVM training is larger than the total number of cells in your dataset. Please adjust svm.train.inds parameter and rerun sc3_process()."
                    )
                )
            }
        }
        # run SVM
        tmp <- prepare_for_svm(ncol(dataset), svm.num.cells, svm.train.inds)
        
        object@consensus$"svm_train_inds" <- tmp$svm.train.inds
        object@consensus$"svm_study_inds" <- tmp$svm.study.inds
    }
    
    # register computing cluster (N-1 CPUs) on a local machine
    if(is.null(n.cores)) {
        n.cores <- parallel::detectCores()
        if(is.null(n.cores)) {
            return("Cannot define a number of available CPU cores that can be used by SC3. Try to set the n.cores parameter in the sc3() function call.")
        }
        # leave one core for the user
        if(n.cores > 1) {
            n.cores <- n.cores - 1
        }
    }
    
    object@consensus$"sc3_n_cores" <- n.cores
    
    if(file.exists(paste0(file.path(find.package("RSelenium"),
                                    "bin/selenium-server-standalone.jar")))) {
        RSelenium::startServer(args=paste("-log", tempfile()), log=FALSE)
        object@consensus$"rselenium" <- TRUE
    } else {
        object@consensus$"rselenium" <- FALSE
    }
    on.exit(stopSeleniumServer())
    
    return(object)
}

#' @export
setMethod("sc3_process", signature(object = "SCESet"),
          function(
              object,
              exprs_values = "counts",
              gene.filter = TRUE,
              gene.filter.fraction = 0.06,
              gene.reads.rare = 2,
              gene.reads.ubiq = 0,
              log.scale = TRUE,
              d.region.min = 0.04,
              d.region.max = 0.07,
              svm.num.cells = NULL,
              svm.train.inds = NULL,
              n.cores = NULL
          ) {
              sc3_process.SCESet(object,
                                exprs_values,
                                gene.filter,
                                gene.filter.fraction,
                                gene.reads.rare,
                                gene.reads.ubiq,
                                log.scale,
                                d.region.min,
                                d.region.max,
                                svm.num.cells,
                                svm.train.inds,
                                n.cores)
          })

#' @export
sc3_calc_dists.SCESet <- function(
    object
) {
    dataset <- object@consensus$sc3_processed_dataset
    if ( is.null(dataset) ) {
        warning(paste0("Please run sc3_process() first!"))
        return(object)
    }
    
    # check whether in the SVM regime
    if(!is.null(object@consensus$svm_train_inds)) {
        dataset <- dataset[ , object@consensus$svm_train_inds]
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- j <- NULL
    
    distances <- c("euclidean", "pearson", "spearman")
    
    message("Calculating distances between the cells...")
    
    if(object@consensus$sc3_n_cores > length(distances)) {
        n.cores <- length(distances)
    } else {
        n.cores <- object@consensus$sc3_n_cores
    }
    
    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    pb <- txtProgressBar(min = 1, max = length(distances), style = 3)
    
    # calculate distances in parallel
    dists <- foreach::foreach(i = distances) %dorng% {
        try({
            setTxtProgressBar(pb, i)
            calculate_distance(dataset, i)
        })
    }
    
    close(pb)
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(dists) <- distances
    
    object@consensus$"sc3_distances" <- dists
    return(object)
}

#' @export
setMethod("sc3_calc_dists", signature(object = "SCESet"),
          function(
              object
          ) {
              sc3_calc_dists.SCESet(object)
          })

#' @export
sc3_calc_transfs.SCESet <- function(
    object
) {
    dists <- object@consensus$sc3_distances
    if ( is.null(dists) ) {
        warning(paste0("Please run sc3_calc_dists() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- j <- NULL
    
    distances <- names(dists)
    transformations <- c("pca", "laplacian")
    
    hash.table <- expand.grid(
        dists = distances,
        transfs = transformations, 
        stringsAsFactors = FALSE
    )
    
    message("Performing transformations...")
    
    if(object@consensus$sc3_n_cores > nrow(hash.table)) {
        n.cores <- nrow(hash.table)
    } else {
        n.cores <- object@consensus$sc3_n_cores
    }
    
    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    pb <- txtProgressBar(min = 1, max = nrow(hash.table), style = 3)
    
    # calculate the 6 distinct transformations in parallel
    transfs <- foreach::foreach(i = 1:nrow(hash.table)) %dopar% {
        try({
            setTxtProgressBar(pb, i)
            transformation(
                get(hash.table[i, 1], dists),
                hash.table[i, 2]
            )
        })
    }
    
    close(pb)
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(transfs) <- paste(hash.table[ , 1], hash.table[ , 2], sep = "_")
    
    object@consensus$"sc3_transformations" <- transfs
    return(object)
}

#' @export
setMethod("sc3_calc_transfs", signature(object = "SCESet"),
          function(
              object
          ) {
              sc3_calc_transfs.SCESet(object)
          })


#' @export
sc3_kmeans.SCESet <- function(
    object,
    ks,
    k.means.iter.max = 1e+09,
    seed = 1
) {
    transfs <- object@consensus$sc3_transformations
    if ( is.null(transfs) ) {
        warning(paste0("Please run sc3_calc_transfs() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    n.dim <- object@consensus$sc3_n_dim

    hash.table <- expand.grid(
        transf = names(transfs),
        ks = ks,
        n.dim = n.dim, 
        stringsAsFactors = FALSE
    )
    
    message("Performing transformations...")
    
    n.cores <- object@consensus$sc3_n_cores

    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    pb <- txtProgressBar(min = 1, max = nrow(hash.table), style = 3)
    
    # calculate the 6 distinct transformations in parallel
    labs <- foreach::foreach(i = 1:nrow(hash.table),
                                .options.RNG = seed) %dopar% {
        try({
            setTxtProgressBar(pb, i)
            transf <- get(hash.table$transf[i], transfs)
            kmeans(
                transf[, 1:hash.table$n.dim[i]],
                hash.table$ks[i],
                iter.max = k.means.iter.max,
                nstart = object@consensus$"sc3_kmeans_nstart"
            )$cluster
        })
    }
    
    close(pb)
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(labs) <- paste(hash.table$transf, hash.table$ks, hash.table$n.dim, sep = "_")
    
    object@consensus$"sc3_kmeans" <- labs
    return(object)
}

#' @export
setMethod("sc3_kmeans", signature(object = "SCESet"),
          function(
              object,
              ks,
              k.means.iter.max = 1e+09,
              seed = 1
          ) {
              sc3_kmeans.SCESet(
                  object,
                  ks,
                  k.means.iter.max,
                  seed)
          })

#' @export
sc3_calc_consens.SCESet <- function(
    object
) {
    k.means <- object@consensus$sc3_kmeans
    if ( is.null(k.means) ) {
        warning(paste0("Please run sc3_kmeans() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    ks <- as.numeric(unique(unlist(lapply(strsplit(names(k.means), "_"), "[[", 3))))
    n.cores <- object@consensus$sc3_n_cores
    
    message("Calculate consensus matrix...")
    
    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    pb <- txtProgressBar(min = 0, max = length(ks), style = 3)

    cons <- foreach::foreach(i = min(ks):max(ks)) %dorng% {
        try({
            setTxtProgressBar(pb, i)
            d <- k.means[grep(paste0("_", i, "_"), names(k.means))]
            d <- unlist(lapply(d, function(x) paste(x, collapse = " ")))

            dat <- consensus_matrix(d)
            tmp <- ED2(dat)
            colnames(tmp) <- as.character(colnames(dat))
            rownames(tmp) <- as.character(colnames(dat))
            diss <- as.dist(as.matrix(as.dist(tmp)))
            hc <- hclust(diss)
            clusts <- get_clusts(hc, i)

            silh <- silhouette(clusts, diss)
            
            labs <- NULL
            for(j in unique(clusts[hc$order])) {
                labs <- rbind(labs, paste(names(clusts[clusts == j]),
                                          collapse = " "))
            }
            
            labs <- as.data.frame(labs)
            colnames(labs) <- "Labels"
            
            list(
                consensus = dat, 
                labels = labs, 
                hc = hc, 
                silhouette = silh
            )
        })
    }
    
    close(pb)
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(cons) <- ks
    
    object@consensus$"sc3_consensus" <- cons
    return(object)
}

#' @export
setMethod("sc3_calc_consens", signature(object = "SCESet"),
          function(
              object
          ) {
              sc3_calc_consens.SCESet(
                  object
              )
          })


#' @export
sc3_calc_biology.SCESet <- function(
    object
) {
    consensus <- object@consensus$sc3_consensus
    if ( is.null(consensus) ) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    ks <- names(consensus)
    n.cores <- object@consensus$sc3_n_cores
    
    message("Calculate consensus matrix...")
    
    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    pb <- txtProgressBar(min = 0, max = length(ks), style = 3)
    
    biol <- foreach::foreach(i = min(ks):max(ks)) %dorng% {
        try({
            hc <- consensus[[as.character(i)]]$hc
            clusts <- get_clusts(hc, i)

            markers <- get_marker_genes(
                object@consensus$sc3_processed_dataset,
                clusts
            )
            
            de.genes <- get_de_genes(
                object@consensus$sc3_processed_dataset,
                clusts
            )
            
            cell.outl <- get_outl_cells(
                object@consensus$sc3_processed_dataset,
                clusts
            )
            
            list(
                markers = markers,
                de.genes = de.genes,
                cell.outl = cell.outl
            )
        })
    }
    
    close(pb)
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(biol) <- ks
    
    object@consensus$"sc3_biology" <- biol
    return(object)
}

#' @export
setMethod("sc3_calc_biology", signature(object = "SCESet"),
          function(
              object
          ) {
              sc3_calc_biology.SCESet(
                  object
              )
          })


#' @export
sc3_run_svm.SCESet <- function(
    object,
    k
) {
    if ( is.null(object@consensus$svm_train_inds) ) {
        warning(paste0("Please rerun sc3_process() defining the training cells!"))
        return(object)
    }
    
    dataset <- object@consensus$sc3_processed_dataset
    hc <- object@consensus$sc3_consensus[[as.character(k)]]$hc
    clusts <- get_clusts(hc, k)

    train.dataset <- dataset[, object@consensus$svm_train_inds]
    colnames(train.dataset) <- clusts
    
    study.labs <-
        support_vector_machines(
            train.dataset,
            dataset[, object@consensus$svm_study_inds],
            "linear"
        )
    
    svm.labs <- c(clusts, study.labs)
    
    ord <- order(c(object@consensus$svm_train_inds, object@consensus$svm_study_inds))
    
    object@consensus$svm_result <- svm.labs[ord]
    return(object)
}

#' @export
setMethod("sc3_run_svm", signature(object = "SCESet"),
          function(
              object,
              k
          ) {
              sc3_run_svm.SCESet(
                  object,
                  k
              )
          })

