#' SC3 main function
#'
#' Run SC3 clustering pipeline on N-1 CPUs and starts the interactive session
#'
#' @param filename either an R matrix / data.frame / data.table object OR a
#' path to your input file containing an expression matrix.
#' @param ks a range of the number of clusters that needs to be tested.
#' k.min is the minimum number of clusters (default is 3). k.max is the maximum
#' number of clusters (default is 7).
#' @param cell.filter defines whether to filter cells that express less than
#' cell.filter.genes genes (lowly expressed cells). By default it is FALSE.
#' Should be used if it is not possible to properly cluster original cells -
#' filtering of lowly expressed cells usually improves clustering.
#' @param cell.filter.genes if cell.filter is used then this parameter defines
#' the minimum number of genes that have to be expressed in each cell
#' (i.e. have more than zero reads). If there are fewer, the cell will be
#' removed from the analysis. The default is 2000.
#' @param gene.filter.fraction the threshold of number of cells to use in the
#' gene filter
#' @param d.region.min the lower boundary of the optimum region of d.
#' The default is 0.04.
#' @param d.region.max the upper boundary of the optimum region of d.
#' The default is 0.07.
#' @param chisq.quantile a treshold used for cell outliers detection.
#' The default is 0.9999.
#' @param interactivity defines whether a browser interactive window should be
#' open after all computation is done. By default it is TRUE. This option can
#' be used to separate clustering calculations from visualisation,
#' e.g. long and time-consuming clustering of really big datasets can be run
#' on a computing cluster and visualisations can be done using a personal
#' laptop afterwards. If interactivity is FALSE then all clustering results
#' will be saved to filename.rds file. To run interactive visulisation with
#' the precomputed clustering results please use
#' sc3_interactive(readRDS(filename.rds)).
#' @param output.directory if interactivity is FALSE, this parameter defines an
#' output directory in which the user has writing rights. The default value is
#' NULL and the output will be written to a temporal directory.
#' @param show.original.labels if cell labels in the dataset are not unique,
#' but represent clusters expected from the experiment, they can be visualised
#' by setting this parameter to TRUE. The default is FALSE.
#' @param svm.num.cells if number of cells in the dataset is greater than this
#' parameter, then an SVM prediction will be used. The default is 1000.
#' @param use.max.cores defines whether to us maximum available cores on the
#' user's machine. Logical, default is TRUE.
#' @param n.cores if use.max.cores is FALSE, this parameter defines the number
#' of cores to be used on the user's machine. Default is 1.
#'
#' @return Opens a browser window with an interactive shine app and visualize
#' all precomputed clusterings.
#'
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom cluster silhouette
#' @importFrom RSelenium startServer
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom utils head write.table setTxtProgressBar txtProgressBar combn
#' @importFrom stats cutree hclust kmeans dist
#'
#' @examples
#' sc3(treutlein, 3:7, interactivity = FALSE, use.max.cores = FALSE)
#'
#' @export
sc3 <- function(filename,
                ks = 3:7,
                cell.filter = FALSE,
                cell.filter.genes = 2000,
                gene.filter.fraction = 0.06,
                d.region.min = 0.04,
                d.region.max = 0.07,
                chisq.quantile = 0.9999,
                interactivity = TRUE,
                output.directory = NULL,
                show.original.labels = FALSE,
                svm.num.cells = 1000,
                use.max.cores = TRUE,
                n.cores = 1) {

    # initial parameters
    set.seed(1)
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    RSelenium::startServer(args=paste("-log", tempfile()), log=FALSE)
    on.exit(stopSeleniumServer())

    # get input data
    dataset <- get_data(filename)

    # remove duplicated genes
    dataset <- dataset[!duplicated(rownames(dataset)), ]

    # cell filter
    if(cell.filter) {
        dataset <- cell_filter(dataset, cell.filter.genes)
    }

    # gene filter
    dataset <- gene_filter(dataset, gene.filter.fraction)

    # log2 transformation
    cat("log2-scaling...\n")
    dataset <- log2(1 + dataset)

    # define the output file basename
    filename <- ifelse(!is.character(filename),
                       deparse(substitute(filename)),
                       basename(filename))

    # define cell names from the input dataset
    cell.names <- c(1:dim(dataset)[2])
    cell.names <- colnames(dataset)

    # prepare for SVM (optional)
    study.cell.names <- NULL
    study.dataset <- data.frame()
    if(dim(dataset)[2] > svm.num.cells) {

        cat("\n")
        cat(paste0("Your dataset contains more than ",
                   svm.num.cells,
                   " cells,therefore clustering wil be performed on a random sample of ",
                   svm.num.cells,
                   " cells, the rest of the cells will be predicted using SVM."))
        cat("\n")
        cat("\n")

        working.sample <- sample(1:dim(dataset)[2], svm.num.cells)
        study.sample <- setdiff(1:dim(dataset)[2], working.sample)
        study.dataset <- dataset[ , study.sample]
        dataset <- dataset[, working.sample]

        study.cell.names <- study.sample
        study.cell.names <- colnames(study.dataset)

        cell.names <- working.sample
        cell.names <- colnames(dataset)
    }

    # define number of cells and region of dimensions
    n.cells <- dim(dataset)[2]
    n.dim <- floor(d.region.min * n.cells) : ceiling(d.region.max * n.cells)

    # for large datasets restrict the region of dimensions to 15
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }

    # create a hash table for running on parallel CPUs
    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              k = c(min(ks) - 1, ks),
                              n.dim = n.dim, stringsAsFactors = FALSE)

    cat("Calculating distance matrices...\n")
    # register computing cluster (N-1 CPUs) on a local machine
    if(use.max.cores) {
        n.cores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)

    # calculate distances in parallel
    dists <- foreach::foreach(i = distances,
                             .options.RNG=1234) %dorng% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances

    # perform kmeans in parallel
    # add a progress bar to be able to see the progress
    pb <- txtProgressBar(min = 1, max = dim(hash.table)[1], style = 3)
    cat("Performing dimensionality reduction and kmeans clusterings...\n")
    labs <- foreach::foreach(i = 1:dim(hash.table)[1],
                            .combine = rbind,
                            .options.RNG=1234) %dorng% {
        try({
            t <- transformation(get(hash.table[i, 1], dists),
                                hash.table[i, 2])[[1]]
            s <- paste(kmeans(t[, 1:hash.table[i, 4]],
                         hash.table[i, 3],
                         iter.max = 1e+09,
                         nstart = 1000)$cluster,
                  collapse = " ")
            setTxtProgressBar(pb, i)
            return(s)
        })
    }
    close(pb)
    res <- cbind(hash.table, labs)
    res$labs <- as.character(res$labs)
    rownames(res) <- NULL

    # perform consensus clustering in parallel
    cat("Computing consensus matrix and labels...\n")
    # first make another hash table for consensus clustering
    all.combinations <- NULL
    for(k in c(min(ks) - 1, ks)) {
        for(i in 1:length(distances)) {
            for(j in 1:length(dimensionality.reductions)) {
                dist.combs <- combn(distances, i)
                dim.red.combs <- combn(dimensionality.reductions, j)
                for(m in 1:dim(dist.combs)[2]) {
                    for(n in 1:dim(dim.red.combs)[2]) {
                        all.combinations <- rbind(
                            all.combinations,
                            cbind(paste(dist.combs[, m], collapse = " "),
                                  paste(dim.red.combs[, n], collapse = " "),
                                  as.numeric(k)))
                    }
                }
            }
        }
    }

    # run consensus clustering in parallel
    cons <- foreach::foreach(i = 1:dim(all.combinations)[1],
                            .options.RNG=1234) %dorng% {
        try({
            d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]] &
                res$k == as.numeric(all.combinations[i, 3]), ]

            dat <- consensus_clustering(d$labs)

            diss <- dist(dat)
            clust <- hclust(diss)
            clusts <- cutree(clust, k = as.numeric(all.combinations[i, 3]))

            silh <- silhouette(clusts, diss)

            labs <- NULL
            for(j in unique(clusts[clust$order])) {
                labs <- rbind(labs, paste(names(clusts[clusts == j]),
                                          collapse = " "))
            }

            labs <- as.data.frame(labs)
            colnames(labs) <- "Labels"

            return(list(dat, labs, clust, silh))
        })
    }

    # stop local cluster
    parallel::stopCluster(cl)

    output.param <- list(filename, distances, dimensionality.reductions,
                         cbind(all.combinations, cons),
                         dataset, study.dataset, svm.num.cells, cell.names,
                         study.cell.names, show.original.labels,
                         chisq.quantile)

    if(interactivity) {
        # start a shiny app in a browser window
        sc3_interactive(output.param)
    } else {
        if(is.null(output.directory)) {
            tmp.dir <- tempdir()
            saveRDS(output.param, paste0(tmp.dir, "/", filename, ".rds"))
            cat(paste0("The output was written to ", tmp.dir, ". Please set the output.directory parameter if you would like to change the output location. Otherwise, please copy the output object ", filename, ".rds", " from the temporal directory before closing this session. The temporal directory will be removed after closing the session."))
        } else {
            saveRDS(output.param, paste0(output.directory, "/", filename, ".rds"))
        }
    }
}
