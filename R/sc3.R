#' SC3 main function
#'
#' Run SC3 clustering pipeline and starts an interactive session in a web browser.
#'
#' @param filename either an R matrix / data.frame object OR a
#' path to your input file containing an input expression matrix. The expression
#' matrix must contain both colnames (cell IDs) and rownames (gene IDs).
#' @param k a range of the number of clusters that needs to be tested.
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
#' @param gene.reads.rare expression value threshold for genes that are expressed in
#' less than gene.filter.fraction*N cells (rare genes)
#' @param gene.reads.ubiq expression value threshold for genes that are expressed in
#' more than (1-fraction)*N cells (ubiquitous genes)
#' @param log.scale defines whether to perform log2 scaling or not. Boolean,
#' default is TRUE.
#' @param d.region.min the lower boundary of the optimum region of d, 
#' default is 0.04.
#' @param d.region.max the upper boundary of the optimum region of d, 
#' default is 0.07.
#' @param k.means.iter.max iter.max parameter used by kmeans() function. Default is 1e+09.
#' @param k.means.nstart nstart parameter used by kmeans() function. Default is 1000.
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
#' @param svm if TRUE then an SVM prediction will be used. The default is FALSE.
#' @param svm.num.cells number of training cells to be used for SVM prediction. 
#' The default is NULL. If the svm parameter is TRUE and svn.num.cells is not provided,
#' then the svm.train.inds parameter is checked.
#' @param svm.train.inds a numeric vector defining indeces of training cells that should be used for SVM training. 
#' The default is NULL. If the svm parameter is TRUE and svn.num.cells and svm.train.inds are not provided,
#' then the defaults of SC3 will be used: if number of cells is more than 5000,
#' then svn.num.cells = 1000, otherwise svn.num.cells = 20 percent of the total number of cells
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
#' @importFrom stats cutree hclust kmeans dist as.dist
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' sc3(treutlein, 3:7, interactivity = FALSE, n.cores = 2)
#'
#' @export
sc3 <- function(dataset,
                k,
                cell.filter = FALSE,
                cell.filter.genes = 2000,
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
  
  # initial parameters
  set.seed(seed)
  distances <- c("euclidean", "pearson", "spearman")
  dimensionality.reductions <- c("PCA", "Spectral")
  if(file.exists(paste0(file.path(find.package("RSelenium"),
                                  "bin/selenium-server-standalone.jar")))) {
    RSelenium::startServer(args=paste("-log", tempfile()), log=FALSE)
    rselenium.installed <- TRUE
  } else {
    rselenium.installed <- FALSE
  }
  on.exit(stopSeleniumServer())

  # get input data
  dataset <- get_data(filename)
  
  # remove duplicated genes
  dataset <- dataset[!duplicated(rownames(dataset)), ]
  
  # cell filter
  if(cell.filter) {
    dataset <- cell_filter(dataset, cell.filter.genes)
    if(ncol(dataset) == 0) {
      cat("All cells were removed after the cell filter! Stopping now...")
      return()
    }
  }
  
  if(ncol(dataset) > 2000) {
      k.means.nstart <- 50
      cat("Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...")
  }
  
    # gene filter
    if(gene.filter) {
        dataset <- gene_filter(dataset, gene.filter.fraction, gene.reads.rare, gene.reads.ubiq)
        if(nrow(dataset) == 0) {
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
  
    # SVM
    
    # prepare for SVM
    study.dataset <- data.frame()
    svm.inds <- NULL
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
                        "You have set both svm.num.cells and svm.train.inds parameters for SVM. Please set only one of them and rerun SC3."
                    )
                )
            }
            if(svm.num.cells >= ncol(dataset) - 1) return(
                message(
                    "Parameter svm.num.cells is larger (or equal) than the number of cells in your dataset. Please make svm.num.cells smaller and rerun SC3."
                )
            )
            if(svm.num.cells < max(ks)) {
                return(
                    message(
                        "Parameter svm.num.cells is smaller than max(ks). Please make svm.num.cells larger and rerun SC3."
                    )
                )
            }
        }
        if(!is.null(svm.train.inds)) {
            if(length(svm.train.inds) < max(ks)) {
                return(
                    message(
                        "The length of svm.train.inds is smaller than the number of clusters k you would like to use for clustering. Please make svm.train.inds parameter longer and rerun SC3."
                    )
                )
            }
            if(max(svm.train.inds) > ncol(dataset) - 1) {
                return(
                    message(
                        "The length of svm.train.inds is larger than the number of cells in your dataset. Please adjust svm.train.inds parameter and rerun SC3."
                    )
                )
            }
        }
        # run SVM
        tmp <- prepare_for_svm(dataset, svm.num.cells, svm.train.inds)
        dataset <- tmp$training.cells
        study.dataset <- tmp$study.cells
        svm.num.cells <- tmp$svm.num.cells
        svm.inds <- tmp$svm.inds
    }
  
  # define number of cells and region of dimensions
  n.cells <- ncol(dataset)
  n.dim <- floor(d.region.min * n.cells) : ceiling(d.region.max * n.cells)
  
  # for large datasets restrict the region of dimensions to 15
  if(length(n.dim) > 15) {
    n.dim <- sample(n.dim, 15)
  }
  
  # create a hash table for running on parallel CPUs
  hash.table <- expand.grid(distan = distances,
                            dim.red = dimensionality.reductions,
                            k = k,
                            n.dim = n.dim, stringsAsFactors = FALSE)
  cat("Calculating distance matrices...\n")
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
  
  cl <- parallel::makeCluster(n.cores, outfile="")
  doParallel::registerDoParallel(cl, cores = n.cores)

  # NULLing the variables to avoid notes in R CMD CHECK
  i <- j <- NULL
  
  # NULLing the variables to avoid notes in R CMD CHECK
  i <- j <- NULL
  
  # calculate distances in parallel
  dists <- foreach::foreach(i = distances) %dorng% {
    try({
      calculate_distance(dataset, i)
    })
  }
  names(dists) <- distances
  
  # add a progress bar to be able to see the progress
  pb <- txtProgressBar(min = 1, max = dim(hash.table)[1], style = 3)
  cat("Performing dimensionality reduction and kmeans clusterings...\n")
  
  # calculate the 6 distinct transformations in parallel
  lis <- foreach::foreach(i = 1:6) %dopar% {
    return(transformation(get(hash.table[i, 1], dists),
                               hash.table[i, 2])[[1]])
  }
  
  # perform kmeans in parallel
  labs <- foreach::foreach(i = 1:dim(hash.table)[1],
                           .combine = rbind,
                           .options.RNG = seed) %dorng% {
                             try({
                               j <- i %% 6
                               if (j == 0) {
                                 j <- 6
                               }
                               t <- lis[[j]]
                               s <- paste(kmeans(t[, 1:hash.table[i, 4]],
                                                 hash.table[i, 3],
                                                 iter.max = k.means.iter.max,
                                                 nstart = k.means.nstart)$cluster,
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
  for(k in k) {
      all.combinations <- rbind(
          all.combinations,
          c(paste(distances, collapse = " "),
            paste(dimensionality.reductions, collapse = " "),
            as.numeric(k))
      )
  }
  # run consensus clustering in parallel
  cons <- foreach::foreach(i = 1:dim(all.combinations)[1]) %dorng% {
    try({
      d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                 res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]] &
                 res$k == as.numeric(all.combinations[i, 3]), ]
      
      dat <- consensus_matrix(d$labs)
      c <- ED2(dat)
      colnames(c) <- as.character(colnames(dat))
      rownames(c) <- as.character(colnames(dat))
      diss <- as.dist(as.matrix(as.dist(c)))
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
  
  output.param <- list(distances = distances,
                       dimensionality.reductions = dimensionality.reductions,
                       cons.table = cbind(all.combinations, cons),
                       dataset = dataset,
                       study.dataset = study.dataset,
                       svm.num.cells = svm.num.cells,
                       svm.inds = svm.inds,
                       show.original.labels = show.original.labels,
                       rselenium.installed = rselenium.installed)
  
    return(output.param)
}
