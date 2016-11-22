#' Run all steps of \code{SC3} in one go
#' 
#' This function is a wrapper that executes all steps of \code{SC3} analysis in one go.
#' 
#' @param object an object of \code{SCESet} class.
#' @param ks a range of the number of clusters \code{k} used for \code{SC3} clustering.
#' Can also be a single integer.
#' @param exprs_values character string 
#' indicating which values should be used
#' as the expression values for \code{SC3} clustering. Valid value is any named element 
#' of the \code{assayData} slot of the \code{SCESet}
#' object. Default is \code{"exprs"}. See \code{\link[scater]{get_exprs}} function of the \code{scater} package for more details.
#' @param gene_filter a boolen variable which defines whether to perform gene 
#' filtering before SC3 clustering.
#' @param pct_dropout_min if \code{gene_filter = TRUE}, then genes with percent of dropouts smaller than 
#' \code{pct_dropout_min} are filtered out before clustering.
#' @param pct_dropout_max if \code{gene_filter = TRUE}, then genes with percent of dropouts larger than 
#' \code{pct_dropout_max} are filtered out before clustering.
#' @param d_region_min defines the minimum number of eigenvectors used for 
#' kmeans clustering as a fraction of the total number of cells. Default is \code{0.04}.
#' See \code{SC3} paper for more details.
#' @param d_region_max defines the maximum number of eigenvectors used for 
#' kmeans clustering as a fraction of the total number of cells. Default is \code{0.07}.
#' See \code{SC3} paper for more details.
#' @param svm_num_cells number of randomly selected training cells to be used 
#' for SVM prediction. The default is \code{NULL}.
#' @param svm_train_inds a numeric vector defining indeces of training cells 
#' that should be used for SVM training. The default is \code{NULL}.
#' @param svm_max define the maximum number of cells below which SVM is not run.
#' @param n_cores defines the number of cores to be used on the user's machine.
#' @param kmeans_nstart nstart parameter passed to \code{\link[stats]{kmeans}} function. Can be set manually. By default it is 
#' \code{1000} for up to \code{2000} cells and \code{50} for more than \code{2000} cells.
#' @param kmeans_iter_max iter.max parameter passed to \code{\link[stats]{kmeans}} 
#' function.
#' @param k_estimator boolean parameter, defines whether to estimate an optimal number of clusters \code{k}.
#' @param biology boolean parameter, defines whether to compute differentially expressed genes, marker 
#' genes and cell outliers.
#' @param rand_seed sets the seed of the random number generator. \code{SC3} is a stochastic
#' method, so setting the \code{rand_seed} to a fixed values can be used for reproducibility
#' purposes.
#' 
#' @name sc3
#' @aliases sc3 sc3,SCESet-method
#' 
#' @return an object of \code{SCESet} class
#' 
#' @export
sc3.SCESet <- function(object, ks = NULL, exprs_values = "exprs", gene_filter = TRUE, pct_dropout_min = 10, 
    pct_dropout_max = 90, d_region_min = 0.04, 
    d_region_max = 0.07, svm_num_cells = NULL, svm_train_inds = NULL, svm_max = 5000, n_cores = NULL, 
    kmeans_nstart = NULL, kmeans_iter_max = 1e+09, k_estimator = FALSE, biology = FALSE, rand_seed = 1) {
    if (is.null(ks)) {
        warning(paste0("Please provide a range of the number of clusters ks to be used by SC3!"))
        return(object)
    }
    object <- sc3_prepare(object, exprs_values, gene_filter, pct_dropout_min, 
        pct_dropout_max, d_region_min, d_region_max, 
        svm_num_cells, svm_train_inds, svm_max, n_cores, kmeans_nstart, kmeans_iter_max, 
        rand_seed)
    if(k_estimator) {
        object <- sc3_estimate_k(object)
    }
    object <- sc3_set_ks(object, ks)
    object <- sc3_calc_dists(object)
    object <- sc3_calc_transfs(object)
    object <- sc3_kmeans(object)
    object <- sc3_calc_consens(object)
    if(biology) {
        object <- sc3_calc_biology(object)
    }
    return(object)
}

#' @rdname sc3
#' @param ... further arguments passed to \code{\link{sc3.SCESet}}
#' @aliases sc3
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3", signature(object = "SCESet"), function(object, ks = NULL, exprs_values = "exprs", 
    gene_filter = TRUE, pct_dropout_min = 10, pct_dropout_max = 90,
    d_region_min = 0.04, d_region_max = 0.07, svm_num_cells = NULL, 
    svm_train_inds = NULL, svm_max = 5000, n_cores = NULL, kmeans_nstart = NULL, kmeans_iter_max = 1e+09, 
    k_estimator = FALSE, biology = FALSE, rand_seed = 1) {
    sc3.SCESet(object, ks, exprs_values, gene_filter, pct_dropout_min, pct_dropout_max, 
        d_region_min, d_region_max, svm_num_cells, svm_train_inds, svm_max, 
        n_cores, kmeans_nstart, kmeans_iter_max, k_estimator, biology, rand_seed)
})

#' Prepare the \code{SCESet} object for \code{SC3} clustering.
#' 
#' This function prepares an object of \code{SCESet} class for \code{SC3} clustering. It
#' creates and populates the following items of the \code{sc3} slot of the \code{SCESet} object:
#' \itemize{
#'   \item \code{exprs_values} - the same as the \code{exprs_values} argument.
#'   \item \code{logged} - a boolen variable which defines whether expression 
#'   values have been log-transformed. If \code{exprs_values != "exprs"} or 
#'   \code{object@logged == FALSE} then it is set to \code{FALSE}. Otherwise it 
#'   is set to \code{TRUE}. Works correctly for all default elements of the
#'   \code{assayData} slot of the \code{SCESet} object. If during the analysis you create your own
#'   element of the \code{assayData} slot, please set the \code{logged} parameter of the 
#'   \code{sc3} slot manually and accordingly after running \code{sc3_prepare}.
#'   \item \code{kmeans_iter_max} - the same as the \code{kmeans_iter_max} argument.
#'   \item \code{kmeans_nstart} - the same as the \code{kmeans_nstart} argument.
#'   \item \code{n_dim} - contains numbers of the number of eigenvectors to be used
#'   in \code{\link[stats]{kmeans}} clustering.
#'   \item \code{rand_seed} - the same as the \code{rand_seed} argument.
#'   \item \code{svm_train_inds} - if SVM is used this item contains indexes of the 
#'   training cells to be used for SC3 clustering and further SVM prediction.
#'   \item \code{svm_study_inds} - if SVM is used this item contains indexes of the
#'    cells to be predicted by SVM.
#'   \item \code{n_cores} - the same as the \code{n_cores} argument.
#'   \item \code{rselenium} - a boolean parameter which defines whether 
#'   \href{https://cran.r-project.org/package=RSelenium}{RSelenium} is properly installed.
#' }
#' 
#' @param object an object of \code{SCESet} class.
#' @param exprs_values character string 
#' indicating which values should be used
#' as the expression values for \code{SC3} clustering. Valid value is any named element 
#' of the \code{assayData} slot of the \code{SCESet}
#' object. Default is \code{"exprs"}. See \code{\link[scater]{get_exprs}} function of the \code{scater} package for more details.
#' @param gene_filter a boolen variable which defines whether to perform gene 
#' filtering before SC3 clustering.
#' @param pct_dropout_min if \code{gene_filter = TRUE}, then genes with percent of dropouts smaller than 
#' \code{pct_dropout_min} are filtered out before clustering.
#' @param pct_dropout_max if \code{gene_filter = TRUE}, then genes with percent of dropouts larger than 
#' \code{pct_dropout_max} are filtered out before clustering.
#' @param d_region_min defines the minimum number of eigenvectors used for 
#' kmeans clustering as a fraction of the total number of cells. Default is \code{0.04}.
#' See \code{SC3} paper for more details.
#' @param d_region_max defines the maximum number of eigenvectors used for 
#' kmeans clustering as a fraction of the total number of cells. Default is \code{0.07}.
#' See \code{SC3} paper for more details.
#' @param svm_num_cells number of randomly selected training cells to be used 
#' for SVM prediction. The default is \code{NULL}.
#' @param svm_train_inds a numeric vector defining indeces of training cells 
#' that should be used for SVM training. The default is \code{NULL}.
#' @param svm_max define the maximum number of cells below which SVM is not run.
#' @param n_cores defines the number of cores to be used on the user's machine.
#' @param kmeans_nstart nstart parameter passed to \code{\link[stats]{kmeans}} function. Default is 
#' \code{1000} for up to \code{2000} cells and \code{50} for more than \code{2000} cells.
#' @param kmeans_iter_max iter.max parameter passed to \code{\link[stats]{kmeans}} 
#' function. Default is \code{1e+09}.
#' @param rand_seed sets the seed of the random number generator. \code{SC3} is a stochastic
#' method, so setting the \code{rand_seed} to a fixed values can be used for reproducibility
#' purposes.
#' 
#' @name sc3_prepare
#' @aliases sc3_prepare sc3_prepare,SCESet-method
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom httr GET
#' @importFrom parallel detectCores
#' @importFrom RSelenium remoteDriver
#' @importFrom scater fData<-
#' @importFrom utils capture.output
#' 
#' @export
sc3_prepare.SCESet <- function(object, exprs_values = "exprs", gene_filter = TRUE, 
    pct_dropout_min = 10, pct_dropout_max = 90,
    d_region_min = 0.04, d_region_max = 0.07, svm_num_cells = NULL, svm_train_inds = NULL, 
    svm_max = 5000, n_cores = NULL, kmeans_nstart = NULL, kmeans_iter_max = 1e+09, rand_seed = 1) {
    
    message("Setting SC3 parameters...")
    
    dataset <- object@assayData[[exprs_values]]
    if (is.null(dataset)) {
        warning(paste0("The object does not contain ", exprs_values, " expression values."))
        return(object)
    }
    
    object@sc3$exprs_values <- exprs_values
    
    object@sc3$logged <- TRUE
    if(exprs_values != "exprs" | !object@logged) {
        object@sc3$logged <- FALSE
    }
    
    # gene filter
    f_data <- object@featureData@data
    f_data$sc3_gene_filter <- TRUE
    if (gene_filter) {
        if(!is.null(f_data$pct_dropout)) {
            f_data$sc3_gene_filter <- f_data$pct_dropout < pct_dropout_max & f_data$pct_dropout > pct_dropout_min
            if (all(!f_data$sc3_gene_filter)) {
                message("All genes were removed after the gene filter! Stopping now...")
                return(object)
            }
        } else {
            warning(paste0("Gene filter can not be calculated, please run calculateQCMetrics() first!"))
            return(object)
        }
    }
    fData(object) <- new("AnnotatedDataFrame", data = f_data)

    object@sc3$kmeans_iter_max <- kmeans_iter_max
    if (is.null(kmeans_nstart)) {
        if (ncol(dataset) > 2000) {
            object@sc3$kmeans_nstart <- 50
            message("Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...")
        } else {
            object@sc3$kmeans_nstart <- 1000
        }
    } else {
        object@sc3$kmeans_nstart <- kmeans_nstart
    }
    
    # define number of cells and region of dimensions
    n_dim <- floor(d_region_min * ncol(dataset)):ceiling(d_region_max * ncol(dataset))
    # for large datasets restrict the region of dimensions to 15
    if (length(n_dim) > 15) {
        n_dim <- sample(n_dim, 15)
    }
    
    # prepare for SVM
    if (!is.null(svm_num_cells) | !is.null(svm_train_inds) | ncol(dataset) > svm_max) {
        # handle all possible errors
        if (!is.null(svm_num_cells)) {
            if (!is.null(svm_train_inds)) {
                return(message("You have set both svm_num_cells and svm_train_inds parameters for SVM training. Please set only one of them and rerun sc3_prepare()."))
            }
            if (svm_num_cells >= ncol(dataset) - 1) 
                return(message("Number of cells used for SVM training is larger (or equal) than the total number of cells in your dataset. Please make svm_num_cells parameter smaller and rerun sc3_prepare()."))
            if (svm_num_cells < 10) {
                return(message("Number of cells used for SVM training is less than 10. Please make sure the number of clusters k is smaller than 10 or increase the number of training cells."))
            }
        }
        if (!is.null(svm_train_inds)) {
            if (length(svm_train_inds) < 10) {
                return(message("Number of cells used for SVM training is less than 10. Please make sure the number of clusters k is smaller than 10 or increase the number of training cells."))
            }
            if (max(svm_train_inds) > ncol(dataset) - 1) {
                return(message("Number of cells used for SVM training is larger than the total number of cells in your dataset. Please adjust svm_train_inds parameter and rerun sc3_prepare()."))
            }
        }
        # run SVM
        tmp <- prepare_for_svm(ncol(dataset), svm_num_cells, svm_train_inds, svm_max)
        
        object@sc3$svm_train_inds <- tmp$svm_train_inds
        object@sc3$svm_study_inds <- tmp$svm_study_inds
        
        # update kmeans_nstart after defining SVM training indeces
        if (is.null(kmeans_nstart)) {
            if (length(tmp$svm_train_inds) <= 2000) {
                object@sc3$kmeans_nstart <- 1000
            }
        } else {
            object@sc3$kmeans_nstart <- kmeans_nstart
        }
        
        # update the region of dimensions
        n_dim <- floor(d_region_min * length(tmp$svm_train_inds)):ceiling(d_region_max * length(tmp$svm_train_inds))
        # for large datasets restrict the region of dimensions to 15
        if (length(n_dim) > 15) {
            n_dim <- sample(n_dim, 15)
        }
    }
    
    object@sc3$n_dim <- n_dim
    
    object@sc3$rand_seed <- rand_seed
    
    # register computing cluster (N-1 CPUs) on a local machine
    if (is.null(n_cores)) {
        n_cores <- parallel::detectCores()
        if (is.null(n_cores)) {
            return("Cannot define a number of available CPU cores that can be used by SC3. Try to set the n_cores parameter in the sc3() function call.")
        }
        # leave one core for the user
        if (n_cores > 1) {
            n_cores <- n_cores - 1
        }
    }
    
    object@sc3$n_cores <- n_cores
    
    # RSelenium is always run on 4444 port
    if (!is.null(tryCatch(httr::GET("http://localhost:4444"), error = function(cond){}))) {
        remDr <- RSelenium::remoteDriver()
        # in the latest version of RSelenium webdriver.gecko.driver has to be set up:
        # http://stackoverflow.com/questions/38751525/firefox-browser-is-not-opening-with-selenium-webbrowser-code
        if(!is.null(tryCatch(utils::capture.output(tmp <- remDr$open()), error = function(cond){}))) {
            object@sc3$rselenium <- TRUE
        } else {
            object@sc3$rselenium <- FALSE
        }
    } else {
        object@sc3$rselenium <- FALSE
    }
    
    return(object)
}

#' @rdname sc3_prepare
#' @aliases sc3_prepare
#' @param ... further arguments passed to \code{\link{sc3_prepare.SCESet}}
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_prepare", signature(object = "SCESet"), function(object, exprs_values = "exprs", 
    gene_filter = TRUE, pct_dropout_min = 10, pct_dropout_max = 90, 
    d_region_min = 0.04, d_region_max = 0.07, svm_num_cells = NULL, 
    svm_train_inds = NULL, svm_max = 5000, n_cores = NULL, kmeans_nstart = NULL, kmeans_iter_max = 1e+09, 
    rand_seed = 1) {
    sc3_prepare.SCESet(object, exprs_values, gene_filter, pct_dropout_min, pct_dropout_max, 
        d_region_min, d_region_max, svm_num_cells, svm_train_inds, 
        svm_max, n_cores, kmeans_nstart, kmeans_iter_max, rand_seed)
})

#' Estimate the optimal k for k-means clustering
#' 
#' Uses Tracy-Widom theory on random matrices to estimate the optimal number of
#' clusters k. Using the function \code{\link{estkTW}} to perform the estimation. 
#' It creates and populates the following items of the `sc3` slot:
#' \itemize{
#'   \item k_estimation - contains the estimated value of `k`.
#' }
#' 
#' @name sc3_estimate_k
#' @aliases sc3_estimate_k sc3_estimate_k,SCESet-method
#' 
#' @param object an object of \code{SCESet} class
#' @return an estimated value of k
#' 
#' @export
sc3_estimate_k.SCESet <- function(object) {
    message("Estimating k...")
    dataset <- get_processed_dataset(object)
    if (is.null(dataset)) {
        warning(paste0("Please run sc3_prepare() first!"))
        return(object)
    }
    res <- estkTW(dataset = dataset)
    object@sc3$k_estimation <- res
    return(object)
}

#' @rdname sc3_estimate_k
#' @aliases sc3_estimate_k
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_estimate_k", signature(object = "SCESet"), function(object) {
    sc3_estimate_k.SCESet(object)
})

#' Sets a range of the number of clusters k used for SC3 clustering.
#' 
#' This function creates and populates the following items of the object@sc3 slot:
#' \itemize{
#'   \item ks - contains a range of the number of clusters k to be used by SC3
#' }
#' 
#' @name sc3_set_ks
#' @aliases sc3_set_ks, sc3_set_ks,SCESet-method
#' 
#' @param object an object of 'SCESet' class
#' @param ks a continuous range of integers - the number of clusters k used for SC3 clustering.
#' Can also be a single integer.
#' 
#' @export
sc3_set_ks.SCESet <- function(object, ks = NULL) {
    if (is.null(ks)) {
        warning(paste0("Please provide a range of the number of clusters ks to be used by SC3!"))
        return(object)
    }
    message("Setting a range of k...")
    object@sc3$ks <- ks
    return(object)
}

#' @rdname sc3_set_ks
#' @aliases sc3_set_ks
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_set_ks", signature(object = "SCESet"), function(object, ks = NULL) {
    sc3_set_ks.SCESet(object, ks)
})

#' Calculate distances between the cells.
#' 
#' This function calculates distances between the cells contained in 
#' the processed_dataset item of the \code{sc3} slot of the \code{SCESet} object. It then
#' creates and populates the following items of the \code{sc3} slot:
#' \itemize{
#'   \item \code{distances} - contains a list of distance matrices corresponding to
#'   Euclidean, Pearson and Spearman distances.
#' }
#' 
#' @name sc3_calc_dists
#' @aliases sc3_calc_dists, sc3_calc_dists,SCESet-method
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @export
sc3_calc_dists.SCESet <- function(object) {
    dataset <- get_processed_dataset(object)
    if (is.null(dataset)) {
        warning(paste0("Please run sc3_prepare() first!"))
        return(object)
    }
    
    # check whether in the SVM regime
    if (!is.null(object@sc3$svm_train_inds)) {
        dataset <- dataset[, object@sc3$svm_train_inds]
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    distances <- c("euclidean", "pearson", "spearman")
    
    message("Calculating distances between the cells...")
    
    if (object@sc3$n_cores > length(distances)) {
        n_cores <- length(distances)
    } else {
        n_cores <- object@sc3$n_cores
    }
    
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    # calculate distances in parallel
    dists <- foreach::foreach(i = distances) %dorng% {
        try({
            calculate_distance(dataset, i)
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(dists) <- distances
    
    object@sc3$distances <- dists
    return(object)
}

#' @rdname sc3_calc_dists
#' @aliases sc3_calc_dists
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_calc_dists", signature(object = "SCESet"), function(object) {
    sc3_calc_dists.SCESet(object)
})

#' Calculate transformations of the distance matrices.
#' 
#' This function transforms all \code{distances} items of the \code{sc3} slot 
#' of the \code{SCESet} object using either principal component analysis (PCA) 
#' or by calculating the eigenvectors of the associated graph Laplacian.
#' The columns of the resulting matrices are then sorted in descending order 
#' by their corresponding eigenvalues. The first \code{d} columns 
#' (where \code{d = max(object@sc3$n_dim)}) of each transformation are then 
#' written to the \code{transformations} item of the \code{sc3} slot.
#' Additionally, this function also removes the previously calculated \code{distances} from
#' the \code{sc3} slot, as they are not needed for further analysis.
#' 
#' @name sc3_calc_transfs
#' @aliases sc3_calc_transfs, sc3_calc_transfs,SCESet-method
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @export
sc3_calc_transfs.SCESet <- function(object) {
    dists <- object@sc3$distances
    if (is.null(dists)) {
        warning(paste0("Please run sc3_calc_dists() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    distances <- names(dists)
    transformations <- c("pca", "laplacian")
    
    n_dim <- object@sc3$n_dim
    
    hash.table <- expand.grid(dists = distances, transfs = transformations, stringsAsFactors = FALSE)
    
    message("Performing transformations and calculating eigenvectors...")
    
    if (object@sc3$n_cores > nrow(hash.table)) {
        n_cores <- nrow(hash.table)
    } else {
        n_cores <- object@sc3$n_cores
    }
    
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    # calculate the 6 distinct transformations in parallel
    transfs <- foreach::foreach(i = 1:nrow(hash.table)) %dopar% {
        try({
            tmp <- transformation(get(hash.table[i, 1], dists), hash.table[i, 2])
            tmp[ , 1:max(n_dim)]
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(transfs) <- paste(hash.table[, 1], hash.table[, 2], sep = "_")
    
    object@sc3$transformations <- transfs
    # remove distances after calculating transformations
    object@sc3$distances <- NULL
    return(object)
}

#' @rdname sc3_calc_transfs
#' @aliases sc3_calc_transfs
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_calc_transfs", signature(object = "SCESet"), function(object) {
    sc3_calc_transfs.SCESet(object)
})

#' \code{kmeans} clustering of cells.
#' 
#' This function performs \code{\link[stats]{kmeans}} clustering of the matrices 
#' contained in the \code{transformations} item of the \code{sc3} slot of the \code{SCESet} object. It then
#' creates and populates the following items of the \code{sc3} slot:
#' \itemize{
#'   \item \code{kmeans} - contains a list of kmeans clusterings.
#' }
#' Additionally, it also removes the previously calculated \code{transformations} from
#' the \code{sc3} slot, as they are not needed for further analysis.
#' 
#' See \code{\link{sc3_prepare}} for the default clustering parameters.
#' 
#' @name sc3_kmeans
#' @aliases sc3_kmeans, sc3_kmeans,SCESet-method
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats kmeans
#' 
#' @export
sc3_kmeans.SCESet <- function(object) {
    transfs <- object@sc3$transformations
    if (is.null(transfs)) {
        warning(paste0("Please run sc3_calc_transfs() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    n_dim <- object@sc3$n_dim
    
    hash.table <- expand.grid(transf = names(transfs), ks = object@sc3$ks, n_dim = n_dim, 
        stringsAsFactors = FALSE)
    
    message("Performing k-means clustering...")
    
    n_cores <- object@sc3$n_cores
    
    kmeans_iter_max <- object@sc3$kmeans_iter_max
    kmeans_nstart <- object@sc3$kmeans_nstart
    rand_seed <- object@sc3$rand_seed
    
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    pb <- utils::txtProgressBar(min = 1, max = nrow(hash.table), style = 3)
    
    # calculate the 6 distinct transformations in parallel
    labs <- foreach::foreach(i = 1:nrow(hash.table), .options.RNG = rand_seed) %dopar% 
        {
            try({
                utils::setTxtProgressBar(pb, i)
                transf <- get(hash.table$transf[i], transfs)
                stats::kmeans(transf[, 1:hash.table$n_dim[i]], hash.table$ks[i], 
                  iter.max = kmeans_iter_max, nstart = kmeans_nstart)$cluster
            })
        }
    
    close(pb)
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(labs) <- paste(hash.table$transf, hash.table$ks, hash.table$n_dim, sep = "_")
    
    object@sc3$kmeans <- labs
    # remove transformations after calculating clusterings
    object@sc3$transformations <- NULL
    return(object)
}

#' @rdname sc3_kmeans
#' @aliases sc3_kmeans
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_kmeans", signature(object = "SCESet"), function(object) {
    sc3_kmeans.SCESet(object)
})

#' Calculate consensus matrix.
#' 
#' This function calculates consensus matrices based on the clustering solutions
#' contained in the \code{kmeans} item of the \code{sc3} slot of the \code{SCESet} object. It then
#' creates and populates the \code{consensus} item of the \code{sc3} slot with 
#' consensus matrices, their hierarchical clusterings in \code{hclust} objects,
#' and Silhouette indeces of the clusters.
#' Additionally, it also removes the previously calculated \code{kmeans} clusterings from
#' the \code{sc3} slot, as they are not needed for further analysis.
#' 
#' @name sc3_calc_consens
#' @aliases sc3_calc_consens, sc3_calc_consens,SCESet-method
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom cluster silhouette
#' @importFrom stats hclust dist as.dist
#' @importFrom scater pData<-
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#' 
#' @export
sc3_calc_consens.SCESet <- function(object) {
    k.means <- object@sc3$kmeans
    if (is.null(k.means)) {
        warning(paste0("Please run sc3_kmeans() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    ks <- as.numeric(unique(unlist(lapply(strsplit(names(k.means), "_"), "[[", 3))))
    
    if (object@sc3$n_cores > length(ks)) {
        n_cores <- length(ks)
    } else {
        n_cores <- object@sc3$n_cores
    }
    
    message("Calculating consensus matrix...")
    
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    cons <- foreach::foreach(i = min(ks):max(ks)) %dorng% {
        try({
            d <- k.means[grep(paste0("_", i, "_"), names(k.means))]
            d <- matrix(unlist(d), nrow=length(d[[1]]))
            dat <- consensus_matrix(d)
            tmp <- ED2(dat)
            colnames(tmp) <- as.character(colnames(dat))
            rownames(tmp) <- as.character(colnames(dat))
            diss <- stats::as.dist(as.matrix(stats::as.dist(tmp)))
            hc <- stats::hclust(diss)
            clusts <- get_clusts(hc, i)
            
            silh <- cluster::silhouette(clusts, diss)
            
            list(consensus = dat, hc = hc, silhouette = silh)
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(cons) <- ks
    
    object@sc3$consensus <- cons
    # remove kmeans results after calculating consensus
    object@sc3$kmeans <- NULL
    
    p_data <- object@phenoData@data
    for(k in ks) {
        hc <- object@sc3$consensus[[as.character(k)]]$hc
        clusts <- get_clusts(hc, k)
        # in case of hybrid SVM approach
        if (!is.null(object@sc3$svm_train_inds)) {
            tmp <- rep(NA, nrow(p_data))
            tmp[object@sc3$svm_train_inds] <- clusts
            clusts <- tmp
        }
        p_data[ , paste0("sc3_", k, "_clusters")] <- clusts
        pData(object) <- new("AnnotatedDataFrame", data = p_data)
    }
    
    return(object)
}

#' @rdname sc3_calc_consens
#' @aliases sc3_calc_consens
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_calc_consens", signature(object = "SCESet"), function(object) {
    sc3_calc_consens.SCESet(object)
})


#' Calculate DE genes, marker genes and cell outliers.
#' 
#' This function calculates DE genes, marker genes and cell outliers based on 
#' the consensus clusterings
#' contained in the consensus item of the object@sc3 slot. It then
#' creates and populates the following items of the object@sc3 slot:
#' \itemize{
#'   \item biology - contains lists of DE genes, marker genes and 
#'   cell outliers data frames.
#' }
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @export
sc3_calc_biology.SCESet <- function(object) {
    consensus <- object@sc3$consensus
    if (is.null(consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    
    message("Computing DE genes, marker genes and cell outliers...")
    
    dataset <- get_processed_dataset(object)
    # check whether in the SVM regime
    if (!is.null(object@sc3$svm_train_inds)) {
        dataset <- dataset[, object@sc3$svm_train_inds]
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    ks <- as.numeric(names(consensus))
    
    if (object@sc3$n_cores > length(ks)) {
        n_cores <- length(ks)
    } else {
        n_cores <- object@sc3$n_cores
    }
    
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    biol <- foreach::foreach(i = min(ks):max(ks)) %dorng% {
        try({
            hc <- consensus[[as.character(i)]]$hc
            clusts <- get_clusts(hc, i)
            
            markers <- get_marker_genes(dataset, clusts)
            
            de.genes <- get_de_genes(dataset, clusts)
            
            cell.outl <- get_outl_cells(dataset, clusts)
            
            list(markers = markers, de.genes = de.genes, cell.outl = cell.outl)
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(biol) <- ks
    
    object@sc3$biology <- biol

    f_data <- object@featureData@data
    p_data <- object@phenoData@data
    for(k in ks) {
        # save DE genes
        f_data[ , paste0("sc3_", k, "_de_padj")] <- NA
        f_data[ , paste0("sc3_", k, "_de_padj")][which(f_data$sc3_gene_filter)] <- biol[[as.character(k)]]$de.genes
        # save marker genes
        f_data[ , paste0("sc3_", k, "_markers_clusts")] <- NA
        f_data[ , paste0("sc3_", k, "_markers_padj")] <- NA
        f_data[ , paste0("sc3_", k, "_markers_auroc")] <- NA
        f_data[ , paste0("sc3_", k, "_markers_clusts")][which(f_data$sc3_gene_filter)] <- biol[[as.character(k)]]$markers[,2]
        f_data[ , paste0("sc3_", k, "_markers_padj")][which(f_data$sc3_gene_filter)] <- biol[[as.character(k)]]$markers[,3]
        f_data[ , paste0("sc3_", k, "_markers_auroc")][which(f_data$sc3_gene_filter)] <- biol[[as.character(k)]]$markers[,1]
        # save cell outliers
        outl <- biol[[as.character(k)]]$cell.outl
        # in case of hybrid SVM approach
        if (!is.null(object@sc3$svm_train_inds)) {
            tmp <- rep(NA, nrow(p_data))
            tmp[object@sc3$svm_train_inds] <- outl
            outl <- tmp
        }
        p_data[ , paste0("sc3_", k, "_log2_outlier_score")] <- log2(outl + 1)
    }
    fData(object) <- new("AnnotatedDataFrame", data = f_data)
    pData(object) <- new("AnnotatedDataFrame", data = p_data)
    return(object)
}

#' @rdname sc3_calc_biology.SCESet
#' @aliases sc3_calc_biology
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_calc_biology", signature(object = "SCESet"), function(object) {
    sc3_calc_biology.SCESet(object)
})


#' Run SVM on training cells
#' 
#' This function performs training of the SVM classifier on the training cells,
#' which indeces are  
#' contained in the svm_train_inds item of the object@sc3 slot. Then it 
#' predicts the labels of the remaining cells using the SVM classifier. Finally it
#' creates and populates the following items of the object@sc3 slot:
#' \itemize{
#'   \item svm_result - contains labels of the cells predicted by the SVM ordered
#'   as the cells in the original dataset.
#' }
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
sc3_run_svm.SCESet <- function(object) {
    if (is.null(object@sc3$svm_train_inds)) {
        warning(paste0("Please rerun sc3_prepare() defining the training cells!"))
        return(object)
    }
    
    dataset <- get_processed_dataset(object)
    
    ks <- names(object@sc3$consensus)
    
    p_data <- object@phenoData@data
    
    for(k in ks) {
        hc <- object@sc3$consensus[[as.character(k)]]$hc
        clusts <- get_clusts(hc, k)
        
        train.dataset <- dataset[, object@sc3$svm_train_inds]
        colnames(train.dataset) <- clusts
        
        study.labs <- support_vector_machines(train.dataset, dataset[, object@sc3$svm_study_inds], 
                                              "linear")
        
        svm.labs <- c(clusts, study.labs)
        
        ord <- order(c(object@sc3$svm_train_inds, object@sc3$svm_study_inds))
        
        p_data[ , paste0("sc3_", k, "_clusters")] <- svm.labs[ord]
        pData(object) <- new("AnnotatedDataFrame", data = p_data)
    }

    return(object)
}

#' @rdname sc3_run_svm.SCESet
#' @aliases sc3_run_svm
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_run_svm", signature(object = "SCESet"), function(object) {
    sc3_run_svm.SCESet(object)
})

#' Summarise SC3 results
#' 
#' This function summarised all SC3 results into a single list and populates 
#' it to the following item of the object@sc3 slot:
#' \itemize{
#'   \item results - contains all SC3 results
#' }
#' 
#' @param object an object of 'SCESet' class
#' @param k the number of clusters k for which the results should be summarised
#' 
#' @importFrom scater pData<-
#' @importFrom methods new
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
sc3_summarise_results.SCESet <- function(object, k) {
    if (is.null(object@sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    p_data <- object@phenoData@data
    if (!paste0("sc3_", k, "_clusters") %in% colnames(p_data)) {
        warning(paste0("There is no clustering result for k = ", k, ", please used a different k!"))
        return(object)
    }
    clusts <- as.data.frame(p_data[, paste0("sc3_", k, "_clusters")])
    colnames(clusts) <- "sc3_clusters"
    rownames(clusts) <- rownames(p_data)
    cell_outliers <- NULL
    de_genes <- NULL
    markers <- NULL
    if (!is.null(object@sc3$biology)) {
        cell_outliers <- object@sc3$biology[[as.character(k)]]$cell.outl
        rownames(cell_outliers) <- rownames(p_data)[!is.na(p_data$sc3_clusters)]
        de_genes <- object@sc3$biology[[as.character(k)]]$de.genes
        markers <- object@sc3$biology[[as.character(k)]]$markers
    }
    res <- list(clusters = clusts, de_genes = de_genes, markers = markers, cell_outliers = cell_outliers)
    object@sc3$results <- res
    return(object)
}

#' @rdname sc3_summarise_results.SCESet
#' @aliases sc3_summarise_results
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_summarise_results", signature(object = "SCESet"), function(object, 
    k) {
    sc3_summarise_results.SCESet(object, k)
})

#' Write SC3 results to Excel file
#' 
#' This function writes SC3 results contained in the object@sc3$results
#' list to an excel file.
#' 
#' @param object an object of 'SCESet' class
#' @param filename name of the excel file, to which the results will be written
#' 
#' @importFrom WriteXLS WriteXLS
#' 
#' @export
sc3_export_results_xls.SCESet <- function(object, filename = "sc3_results.xls") {
    if (is.null(object@sc3$results)) {
        stop(paste0("Please run sc3_summarise_results() first!"))
    }
    WriteXLS(object@sc3$results, ExcelFileName = filename, SheetNames = names(object@sc3$results), 
        row.names = TRUE)
}

#' @rdname sc3_export_results_xls.SCESet
#' @aliases sc3_export_results_xls
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_export_results_xls", signature(object = "SCESet"), function(object, 
    filename = "sc3_results.xls") {
    sc3_export_results_xls.SCESet(object, filename)
})
