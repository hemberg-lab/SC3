#' Single cell RNA-Seq data extracted from a publication by Treutlein et al.
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52583}
#'
#' Columns represent cells, rows represent genes expression values. Colnames
#' respresent indexes of cell clusters (known information based on the
#' experimental protocol). There are 80 cells and 5 clusters in this dataset.
#'
"treutlein"

#' Gene filter
#'
#' The gene filter removes genes/transcripts that are either expressed 
#' (expression value is more than 2) in less than X% of cells 
#' (rare genes/transcripts) or expressed (expression value is more than 0) 
#' in at least (100-X)% of cells (ubiquitous genes/transcripts).
#'
#' @param data input expression matrix
#' @param fraction fraction of cells (X/100).
#' @param reads.rare expression value threshold for genes that are expressed in
#' less than fraction*N cells (rare genes)
#' @param reads.ubiq expression value threshold for genes that are expressed in
#' more than (1-fraction)*N cells (ubiquitous genes)
#' @return filtered expression matrix some genes were removed.
gene_filter <- function(data, fraction = 0.06, reads.rare = 2, reads.ubiq = 0) {
    message("Gene filtering...")
    frac.cells <- ceiling(fraction * ncol(data))
    d <- data[rowSums(data > reads.rare) >= frac.cells & rowSums(data > reads.ubiq) <= 
        ncol(data) - frac.cells, ]
    d <- unique(d)
    return(d)
}

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
        L <- norm_laplacian(exp(-dists/max(dists)), 0)
        l <- eigen(L)
        # sort eigenvectors by their eigenvalues
        return(l$vectors[, order(l$values)])
    } else if (method == "laplacian_reg") {
        L <- norm_laplacian(exp(-dists/max(dists)), 1000)
        l <- eigen(L)
        # sort eigenvectors by their eigenvalues in increasing order
        return(l$vectors[, order(l$values)])
    } else if (method == "mds") {
        t <- cmdscale(dists, k = ncol(dists) - 1)
        return(t[[1]])
    }
}

#' Graph Laplacian calculation
#'
#' Calculate graph Laplacian of a distance matrix
#'
#' @param x adjacency/distance matrix
#' @param tau regularization term
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#' 
#' @return graph Laplacian of the adjacency/distance matrix
norm_laplacian <- function(x, tau) {
    D <- diag(colSums(x)^(-0.5))
    dim <- nrow(x)
    return(diag(dim(D)[1]) - mult(D, x, dim))
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
#' @param clusts a list of clustering solutions. Labels in each solutions 
#' should be separated by a space.
#' @return consensus matrix
#'
#' @importFrom stats dist
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#' @export
consensus_matrix <- function(clusts) {
    n.cells <- length(unlist(strsplit(clusts[1], " ")))
    res <- matrix(0, nrow = n.cells, ncol = n.cells)
    j <- length(clusts)
    res <- consmx(clusts, res, j)
    res <- -res/j
    colnames(res) <- as.character(c(1:n.cells))
    rownames(res) <- as.character(c(1:n.cells))
    return(res)
}

#' Run support vector machines (SVM) prediction
#'
#' Train an SVM classifier on training cells and then
#' classify study cells using the classifier.
#'
#' @param train expression matrix with training cells
#' @param study expression matrix with study cells
#' @param kern kernel to be used with SVM
#' @return classification of study cells
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
#' Defines train and study cell indeces based on the svm.num.cells and
#' svm.train.inds input parameters
#'
#' @param N number of cells in the input dataset
#' @param svm.num.cells number of random cells to be used for training
#' @param svm.train.inds indeces of cells to be used for training
#' @return A list of indeces of the train and the study cells
prepare_for_svm <- function(N, svm.num.cells = NULL, svm.train.inds = NULL) {
    
    if (!is.null(svm.num.cells)) {
        message("Defining training cells for SVM using svm.num.cells parameter...")
        train.inds <- sample(1:N, svm.num.cells)
        study.inds <- setdiff(1:N, train.inds)
    }
    
    if (!is.null(svm.train.inds)) {
        message("Defining training cells for SVM using svm.train.inds parameter...")
        train.inds <- svm.train.inds
        study.inds <- setdiff(1:N, svm.train.inds)
    }
    
    if (is.null(svm.num.cells) & is.null(svm.train.inds)) {
        message("Defining training cells for SVM using 5000 random cells...")
        train.inds <- sample(1:N, 5000)
        study.inds <- setdiff(1:N, train.inds)
    }
    
    return(list(svm.train.inds = train.inds, svm.study.inds = study.inds))
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
    
    message("Computing eigenvalues...")
    
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
