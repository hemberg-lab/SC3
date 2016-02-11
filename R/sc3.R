#' SC3 main function
#'
#' Run SC3 clustering pipeline and starts an interactive session in a web browser.
#'
#' @param filename either an R matrix / data.frame object OR a
#' path to your input file containing an input expression matrix. The expression
#' matrix must contain both colnames (cell IDs) and rownames (gene IDs).
#' @param ks a range of the number of clusters that needs to be tested.
#' k.min is the minimum number of clusters (default is 3). k.max is the maximum
#' number of clusters (default is 7).
#' @param cell.filter defines whether to filter cells that express less than
#' cell.filter.genes genes (lowly expressed cells). By default it is FALSE.
#' The cell filter should be used if the quality of data is low, i.e. if one
#' suspects that some of the cells may be technical outliers with poor coverage.
#' Filtering of lowly expressed cells usually improves clustering.
#' @param cell.filter.genes if cell.filter is used then this parameter defines
#' the minimum number of genes that have to be expressed in each cell
#' (expression value > 1e-2). If there are fewer, the cell will be
#' removed from the analysis. The default is 2000.
#' @param gene.filter defines whether to perform gene filtering or not. Boolean,
#' default is TRUE.
#' @param gene.filter.fraction fraction of cells (1 - X/100), default is 0.06.
#' The gene filter removes genes that are either expressed or absent
#' (expression value is less than 2) in at least X % of cells.
#' The motivation for the gene filter is that ubiquitous and rare genes most
#' often are not informative for the clustering.
#' @param log.scale defines whether to perform log2 scaling or not. Boolean,
#' default is TRUE.
#' @param d.region.min the lower boundary of the optimum region of d, 
#' default is 0.04.
#' @param d.region.max the upper boundary of the optimum region of d, 
#' default is 0.07.
#' @param chisq.quantile a threshold used for cell outliers detection, 
#' default is 0.9999.
#' @param auroc.threshold a threshold used for marker genes detection,
#' default is 0.85
#' @param interactivity defines whether a browser interactive window should be
#' open after all computation is done. By default it is TRUE. This option can
#' be used to separate clustering calculations from visualisation,
#' e.g. long and time-consuming clustering of really big datasets can be run
#' on a farm cluster and visualisations can be done using a personal
#' laptop afterwards. If interactivity is FALSE then all clustering results
#' will be saved as "sc3.interactive.arg" list. To run interactive visulisation with
#' the precomputed clustering results please use
#' sc3_interactive(sc3.interactive.arg).
#' @param show.original.labels if cell labels in the dataset are not unique,
#' but represent clusters expected from the experiment, they can be visualised
#' by setting this parameter to TRUE. The default is FALSE.
#' @param svm.num.cells if number of cells in the dataset is greater than this
#' parameter, then an SVM prediction will be used. The default is 1000.
#' @param use.max.cores defines whether to use maximum available cores on the
#' user's machine. Logical, default is TRUE.
#' @param n.cores defines the number
#' of cores to be used on the user's machine. Default is NA.
#' @param seed sets seed for the random number generator, default is 1.
#' Can be used to check the stability of clustering results: if the results are 
#' the same after changing the seed several time, then the clustering solution 
#' is stable.
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
#' sc3(treutlein, 3:7, interactivity = FALSE, n.cores = 1)
#'
#' @export
sc3 <- function(filename,
                ks = 3:7,
                cell.filter = FALSE,
                cell.filter.genes = 2000,
                gene.filter = TRUE,
                gene.filter.fraction = 0.06,
                log.scale = TRUE,
                d.region.min = 0.04,
                d.region.max = 0.07,
                chisq.quantile = 0.9999,
                auroc.threshold = 0.85,
                interactivity = TRUE,
                show.original.labels = FALSE,
                svm.num.cells = 1000,
                use.max.cores = TRUE,
                n.cores = NA,
                seed = 1) {

    # initial parameters
    set.seed(seed)
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
        if(dim(dataset)[2] == 0) {
            cat("All cells were removed after the cell filter! Stopping now...")
            return()
        }
    }

    # gene filter
    if(gene.filter) {
        dataset <- gene_filter(dataset, gene.filter.fraction)
        if(dim(dataset)[1] == 0) {
            cat("All genes were removed after the gene filter! Stopping now...")
            return()
        }
    }

    # log2 transformation
    if(log.scale) {
        cat("log2-scaling...\n")
        dataset <- log2(1 + dataset)
    }

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
    if(is.na(n.cores)) {
        if(use.max.cores) {
            n.cores <- parallel::detectCores() - 1
        } else {
            return("Please provide a number of CPU cores (n.cores) that can be used by SC3.")
        }
    }
    
    cl <- parallel::makeCluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, cores = n.cores)

    # calculate distances in parallel
    dists <- foreach::foreach(i = distances) %dorng% {
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
                            .options.RNG = seed) %dorng% {
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
    cons <- foreach::foreach(i = 1:dim(all.combinations)[1]) %dorng% {
        try({
            d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]] &
                res$k == as.numeric(all.combinations[i, 3]), ]

            dat <- consensus_matrix(d$labs)

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

    output.param <- list(filename = filename,
                         distances = distances,
                         dimensionality.reductions = dimensionality.reductions,
                         cons.table = cbind(all.combinations, cons),
                         dataset = dataset,
                         study.dataset = study.dataset,
                         svm.num.cells = svm.num.cells,
                         cell.names = cell.names,
                         study.cell.names = study.cell.names,
                         show.original.labels = show.original.labels,
                         chisq.quantile = chisq.quantile,
                         auroc.threshold = auroc.threshold)

    if(interactivity) {
        # start a shiny app in a browser window
        sc3_interactive(output.param)
    } else {
        sc3.interactive.arg <- list()
        sc3.interactive.arg <<- output.param
    }
}
