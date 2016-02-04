#' Import expression matrix
#'
#' Import an input expression matrix.
#'
#' @param name name of an R object or a text file
#' @return expression matrix
#'
#' @importFrom utils read.table read.csv
get_data <- function(name) {
    if(!is.character(name)) {
        return(name)
    } else {
        if(!grepl("csv", name)) {
            return(as.matrix(read.table(name, check.names = FALSE)))
        } else if(grepl("csv", name)) {
            return(as.matrix(read.csv(name, check.names = FALSE)))
        }
    }
}

#' Cell filter
#'
#' The cell filter should be used if the quality of data is low, i.e. if one
#' suspects that some of the cells may be technical outliers with poor coverage.
#' The cell filter removes cells containing fewer than cell.filter.genes.
#'
#' @param data expression matrix
#' @param cell.filter.genes minimum number of genes that must be expressed in
#' each cell, default is 2,000.
#' @return Filtered expression matrix
cell_filter <- function(data, cell.filter.genes) {
    cat("Cell filtering...\n")
    data <- data[ , colSums(data > 1e-2) > cell.filter.genes]
    return(data)
}

#' Gene filter
#'
#' The gene filter removes genes that are either expressed or absent
#' (expression value is less than 2) in at least X % of cells.
#' The motivation for the gene filter is that ubiquitous and rare genes most
#' often are not informative for the clustering.
#'
#' @param data expression matrix
#' @param fraction fraction of cells (1 - X/100), default is 0.06.
#' @return filtered expression matrix some genes were removed.
gene_filter <- function(data, fraction) {
    cat("Gene filtering...\n")
    filter.params <- filter_params(data, fraction)
    min.cells <- filter.params$min.cells
    max.cells <- filter.params$max.cells
    min.reads <- filter.params$min.reads
    d <- data[rowSums(data > min.reads) >= min.cells &
              rowSums(data > 0) <= dim(data)[2] - max.cells, ]
    d <- unique(d)
    return(d)
}

filter_params <- function(dataset, fraction) {
    n.cells <- dim(dataset)[2]

    min.cells <- ceiling(fraction*n.cells)
    max.cells <- ceiling(fraction*n.cells)
    min.reads <- 2

    return(list(min.cells = min.cells, max.cells = max.cells,
                min.reads = min.reads))
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
        as.matrix(dist(t(data), method = method))
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
    if (method == "pca") {
        t <- prcomp(dists, center = TRUE, scale. = TRUE)
        list(t$rotation, t$sdev)
    } else if (method == "spectral") {
        L <- norm_laplacian(exp(-dists/max(dists)), 0)
        # sort eigenvectors by their eigenvalues in increasing order
        list(eigen(L)$vectors[, order(eigen(L)$values)],
             eigen(L)$values[order(eigen(L)$values)])
    } else if (method == "spectral_reg") {
        L <- norm_laplacian(exp(-dists/max(dists)), 1000)
        # sort eigenvectors by their eigenvalues in increasing order
        list(eigen(L)$vectors[, order(eigen(L)$values)],
             eigen(L)$values[order(eigen(L)$values)])
    } else if (method == "mds") {
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
    x <- x + tau * matrix(1, dim(x)[1], dim(x)[2])
    D <- diag(colSums(x))
    D1 <- D^(-0.5)
    D1[D1 == Inf] <- 0
    return(diag(dim(D)[1]) - D1 %*% x %*% D1)
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
    for (i in 1:length(clusts)) {
        t <- clusts[i]
        t <- as.numeric(unlist(strsplit(t, " ")))
        t <- as.matrix(dist(t))
        t[t != 0] <- -1
        t[t == 0] <- 1
        t[t == -1] <- 0
        res <- res + t
    }
    res <- res/i
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
