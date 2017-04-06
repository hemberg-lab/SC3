#' Single cell RNA-Seq data extracted from a publication by Treutlein et al.
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52583}
#'
#' Columns represent cells, rows represent genes expression values. Colnames
#' respresent indexes of cell clusters (known information based on the
#' experimental protocol). There are 80 cells and 5 clusters in this dataset.
#'
"treutlein"

#' Calculate a distance matrix
#'
#' Distance between the cells, i.e. columns, in the input expression matrix are
#' calculated using the Euclidean, Pearson and Spearman metrics to construct
#' distance matrices.
#'
#' @param data expression matrix
#' @param method one of the distance metrics: 'spearman', 'pearson', 'euclidean'
#' @return distance matrix
#'
#' @importFrom stats cor dist
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#'
calculate_distance <- function(data, method) {
    return(if (method == "spearman") {
        as.matrix(1 - cor(data, method = "spearman"))
    } else if (method == "pearson") {
        as.matrix(1 - cor(data, method = "pearson"))
    } else {
        ED2(data)
    })
}

#' Distance matrix transformation
#'
#' All distance matrices are transformed using either principal component 
#' analysis (PCA) or by calculating the 
#' eigenvectors of the graph Laplacian (Spectral). 
#' The columns of the resulting matrices are then sorted in 
#' descending order by their corresponding eigenvalues.
#'
#' @param dists distance matrix
#' @param method transformation method: either 'pca' or
#' 'laplacian'
#' @return transformed distance matrix
#'
#' @importFrom stats prcomp cmdscale
#'
transformation <- function(dists, method) {
    if (method == "pca") {
        t <- prcomp(dists, center = TRUE, scale. = TRUE)
        return(t$rotation)
    } else if (method == "laplacian") {
        L <- norm_laplacian(dists)
        l <- eigen(L)
        # sort eigenvectors by their eigenvalues
        return(l$vectors[, order(l$values)])
    }
}

#' Calculate consensus matrix
#'
#' Consensus matrix is calculated using the Cluster-based Similarity 
#' Partitioning Algorithm (CSPA). For each clustering solution a binary 
#' similarity matrix is constructed from the corresponding cell labels: 
#' if two cells belong to the same cluster, their similarity is 1, otherwise 
#' the similarity is 0. A consensus matrix is calculated by averaging all 
#' similarity matrices.
#'
#' @param clusts a matrix containing clustering solutions in columns
#' @return consensus matrix
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#' @export
consensus_matrix <- function(clusts) {
    res <- consmx(clusts)
    colnames(res) <- as.character(c(1:nrow(clusts)))
    rownames(res) <- as.character(c(1:nrow(clusts)))
    return(res)
}

#' Run support vector machines (\code{SVM}) prediction
#'
#' Train an \code{SVM} classifier on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training dataset with colnames, corresponding to training labels
#' @param study study dataset
#' @param kern kernel to be used with SVM
#' @return classification of the study dataset
#'
#' @importFrom e1071 svm
#' @importFrom stats predict
support_vector_machines <- function(train, study, kern) {
    train <- t(train)
    labs <- factor(rownames(train))
    rownames(train) <- NULL
    model <- tryCatch(svm(train, labs, kernel = kern), error = function(cond) return(NA))
    pred <- predict(model, t(study))
    return(pred = pred)
}

#' A helper function for the SVM analysis
#'
#' Defines train and study cell indeces based on the svm_num_cells and
#' svm_train_inds input parameters
#'
#' @param N number of cells in the input dataset
#' @param svm_num_cells number of random cells to be used for training
#' @param svm_train_inds indeces of cells to be used for training
#' @param svm_max define the maximum number of cells below which SVM is not run
#' @return A list of indeces of the train and the study cells
prepare_for_svm <- function(N, svm_num_cells = NULL, svm_train_inds = NULL, svm_max) {
    
    if (!is.null(svm_num_cells)) {
        message("Defining training cells for SVM using svm_num_cells parameter...")
        train_inds <- sample(1:N, svm_num_cells)
        study_inds <- setdiff(1:N, train_inds)
    }
    
    if (!is.null(svm_train_inds)) {
        message("Defining training cells for SVM using svm_train_inds parameter...")
        train_inds <- svm_train_inds
        study_inds <- setdiff(1:N, svm_train_inds)
    }
    
    if (is.null(svm_num_cells) & is.null(svm_train_inds)) {
        message(paste0("Defining training cells for SVM using ", svm_max, " random cells..."))
        train_inds <- sample(1:N, svm_max)
        study_inds <- setdiff(1:N, train_inds)
    }
    
    return(list(svm_train_inds = train_inds, svm_study_inds = study_inds))
}

#' Estimate the optimal k for k-means clustering
#' 
#' The function finds the eigenvalues of the sample covariance matrix. 
#' It will then return the number of significant eigenvalues according to 
#' the Tracy-Widom test.
#' 
#' @param dataset processed input expression matrix.
#' @return an estimated number of clusters k
estkTW <- function(dataset) {
    
    p <- ncol(dataset)
    n <- nrow(dataset)
    
    # compute Tracy-Widom bound
    x <- scale(dataset)
    muTW <- (sqrt(n - 1) + sqrt(p))^2
    sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
    sigmaHatNaive <- tmult(x)  # x left-multiplied by its transpose
    bd <- 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
    
    # compute eigenvalues and return the amount which falls above the bound
    evals <- eigen(sigmaHatNaive, symmetric = TRUE, only.values = TRUE)$values
    k <- 0
    for (i in 1:length(evals)) {
        if (evals[i] > bd) {
            k <- k + 1
        }
    }
    return(k)
}

make_col_ann_for_heatmaps <- function(object, show_pdata) {
    if (any(!show_pdata %in% colnames(object@phenoData@data))) {
        show_pdata_excl <- show_pdata[!show_pdata %in% colnames(object@phenoData@data)]
        show_pdata <- show_pdata[show_pdata %in% colnames(object@phenoData@data)]
        message(paste0("Provided columns '", paste(show_pdata_excl, collapse = "', '"), "' do not exist in the phenoData table!"))
        if (length(show_pdata) == 0) {
            return(NULL)
        }
    }
    ann <- NULL
    if (is.null(object@sc3$svm_train_inds)) {
        ann <- object@phenoData@data[, colnames(object@phenoData@data) %in% show_pdata]
    } else {
        ann <- object@phenoData@data[object@sc3$svm_train_inds, colnames(object@phenoData@data) %in% 
            show_pdata]
    }
    # remove columns with 1 value only
    if (length(show_pdata) > 1) {
        keep <- unlist(lapply(ann, function(x) {
            length(unique(x))
        })) > 1
        if (!all(keep)) {
            message(paste0("Columns '", paste(names(keep)[!keep], collapse = "', '"), "' were excluded from annotation since they contained only a single value."))
        }
        ann <- ann[, names(keep)[keep]]
        if (ncol(ann) == 0) {
            ann <- NULL
        } else {
            ann <- as.data.frame(lapply(ann, function(x) {
                if (nlevels(as.factor(x)) > 9) 
                  x else as.factor(x)
            }))
            # convert outlier scores back to numeric
            for (i in grep("_log2_outlier_score", colnames(ann))) {
                if (class(ann[, i]) == "factor") {
                  ann[, i] <- as.numeric(levels(ann[, i]))[ann[, i]]
                }
            }
        }
    } else {
        if (length(unique(ann)) > 1) {
            ann <- as.data.frame(ann)
            colnames(ann) <- show_pdata
            if (!grepl("_log2_outlier_score", show_pdata)) {
                ann <- as.data.frame(lapply(ann, function(x) {
                  if (nlevels(as.factor(x)) > 9) 
                    return(x) else return(as.factor(x))
                }))
            }
        } else {
            message(paste0("Column '", show_pdata, "' was excluded from annotation since they contained only a single value."))
            ann <- NULL
        }
    }
    return(ann)
}

#' Get processed dataset used by SC3 from the default scater slots
#' 
#' Takes data from the 'exprs' slot and applies the gene filter
#' 
#' @param object an object of 'SCESet' class
#' 
#' @importFrom scater get_exprs
#'
#' @export
get_processed_dataset <- function(object) {
    dataset <- scater::get_exprs(object, "exprs")
    if (!is.null(object@featureData@data$sc3_gene_filter)) {
        dataset <- dataset[object@featureData@data$sc3_gene_filter, ]
    }
    return(dataset)
}
