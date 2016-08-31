
#' Gene filter
#'
#' The gene filter removes genes that are either expressed or absent
#' (expression value is less than 2) in at least X % of cells.
#' The motivation for the gene filter is that ubiquitous and rare genes most
#' often are not informative for the clustering.
#'
#' @param data expression matrix
#' @param fraction fraction of cells (1 - X/100), default is 0.06.
#' @param reads.rare expression value threshold for genes that are expressed in
#' less than fraction*N cells (rare genes)
#' @param reads.ubiq expression value threshold for genes that are expressed in
#' more than (1-fraction)*N cells (ubiquitous genes)
#' @return filtered expression matrix some genes were removed.
gene_filter <- function(data, fraction = 0, reads.rare = 0, reads.ubiq = 0) {
    cat("Gene filtering...\n")
    frac.cells <- ceiling(fraction*ncol(data))
    d <- data[rowSums(data > reads.rare) >= frac.cells &
                  rowSums(data > reads.ubiq) <= ncol(data) - frac.cells, ]
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
#' @param method one of the distance metrics: "spearman", "pearson", "euclidean",
#' "maximum", "manhattan", "canberra", "binary" or "minkowski"
#' @return distance matrix
#'
#' @importFrom stats cor dist
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
#' analysis (PCA), multidimensional scaling (MDS) or by calculating the 
#' eigenvectors of the graph Laplacian (Spectral). 
#' The columns of the resulting matrices are then sorted in 
#' descending order by their corresponding eigenvalues.
#'
#' @param dists distance matrix
#' @param method transformation method: either "pca", "mds", "spectral" or
#' "spectral_reg", where "spectral_reg" calculates graph Laplacian with
#' regularization (tau = 1000)
#' @return transformed distance matrix
#'
#' @importFrom stats prcomp cmdscale
#'
transformation <- function(dists, method) {
  if (method == "PCA") {
    t <- prcomp(dists, center = TRUE, scale. = TRUE)
    list(t$rotation, t$sdev)
  } else if (method == "Spectral") {
    L <- norm_laplacian(exp(-dists/max(dists)), 0)
    # sort eigenvectors by their eigenvalues in increasing order
    l <- eigen(L)
    list(l$vectors[, order(l$values)],
         l$values[order(l$values)])
  } else if (method == "spectral_reg") {
    L <- norm_laplacian(exp(-dists/max(dists)), 1000)
    # sort eigenvectors by their eigenvalues in increasing order
    list(eigen(L)$vectors[, order(eigen(L)$values)],
         eigen(L)$values[order(eigen(L)$values)])
  } else if (method == "MDS") {
    t <- cmdscale(dists, k = ncol(dists) - 1)
    list(t)
  }
}

#' Graph Laplacian calculation
#'
#' Calculate graph Laplacian of a distance matrix
#'
#' @param x adjacency/distance matrix
#' @param tau regularization term
#' @return graph Laplacian of the adjacency/distance matrix
norm_laplacian <- function(x, tau) {
  D <- diag(colSums(x)^(-0.5))
  dim <- nrow(x)
  return(diag(dim(D)[1]) - mult(D, x, dim))
}

#' Calculate consensus matrix
#'
#' Consensus matrix is calculated using the Cluster-based Similarity 
#' Partitioning Algorithm (CSPA). For each clustering result a binary 
#' similarity matrix is constructed from the corresponding cell labels: 
#' if two cells belong to the same cluster, their similarity is 1, otherwise 
#' the similarity is 0. A consensus matrix is calculated by averaging all 
#' similarity matrices.
#'
#' @param clusts a list clustering labels (separated by a space)
#' @return consensus matrix
#'
#' @importFrom stats dist
#'
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
  model <- tryCatch(svm(train, labs, kernel = kern),
                    error = function(cond) return(NA))
  pred <- predict(model, t(study))
  return(pred = pred)
}

prepare_for_svm <- function(dataset, svm.num.cells, svm.train.inds) {
    
    svm.inds <- NULL
    
    if(!is.null(svm.num.cells)) {
        cat("Defining training cells for SVM using svm.num.cells parameter...\n")
        svm.train.inds <- sample(1:ncol(dataset), svm.num.cells)
        svm.study.inds <- setdiff(1:ncol(dataset), svm.train.inds)
        study.dataset <- dataset[ , svm.study.inds]
        dataset <- dataset[, svm.train.inds]
        svm.inds <- c(svm.train.inds, svm.study.inds)
    }
    
    if(is.null(svm.inds) & !is.null(svm.train.inds)) {
        cat("Defining training cells for SVM using svm.train.inds parameter...\n")
        svm.num.cells <- length(svm.train.inds)
        svm.study.inds <- setdiff(1:ncol(dataset), svm.train.inds)
        study.dataset <- dataset[ , svm.study.inds]
        dataset <- dataset[, svm.train.inds]
        svm.inds <- c(svm.train.inds, svm.study.inds)
    }
    
    if(is.null(svm.inds) & ncol(dataset) > 5000) {
        cat("Defining training cells for SVM using 5000 random cells...\n")
        svm.num.cells <- 5000
        svm.train.inds <- sample(1:ncol(dataset), svm.num.cells)
        svm.study.inds <- setdiff(1:ncol(dataset), svm.train.inds)
        study.dataset <- dataset[ , svm.study.inds]
        dataset <- dataset[, svm.train.inds]
        svm.inds <- c(svm.train.inds, svm.study.inds)
    }
    
    return(
        list(
            training.cells = dataset, 
            study.cells = study.dataset,
            svm.num.cells = svm.num.cells,
            svm.inds = svm.inds
        )
    )
}
