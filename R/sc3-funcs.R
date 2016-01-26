#' Import expression matrix
#'
#' Import either treutlein dataset or a dataset defined by a user
#'
#' @param name either treutlein or a name of a text file
#' @return Expression matrix
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
#' More than cell.filter.genes have to be expressed in each cell
#'
#' @param data Input expression matrix
#' @param cell.filter.genes number of genes that have to expressed in each cell
#' @return Filtered expression matrix
cell_filter <- function(data, cell.filter.genes) {
    cat("Cell filtering...\n")
    data <- data[ , colSums(data > 1e-2) > cell.filter.genes]
    if(dim(data)[2] == 0) {
        cat(paste0("Your dataset did not pass cell filter (more than ",
                   cell.filter.genes,
                   " genes have to be expressed in each cell)!
                   Stopping now..."))
        return()
    } else {
        return(data)
    }
}

#' Gene filter
#'
#' Filter genes that would not contribute to clustering, because they are either
#' expressed or not expressed in almost all cells.
#'
#' @param data Input expression matrix
#' @param fraction threshold defining a fraction cells for each gene
#' @return Filtered expression matrix
gene_filter <- function(data, fraction) {
    cat("Gene filtering...\n")
    filter.params <- filter_params(data, fraction)
    min.cells <- filter.params$min.cells
    max.cells <- filter.params$max.cells
    min.reads <- filter.params$min.reads
    data <- gene_filter1(data, min.cells, max.cells, min.reads)

    if(dim(data)[1] == 0) {
        cat("All genes were removed after the gene filter! Stopping now...")
        return()
    } else {
        return(data)
    }
}

filter_params <- function(dataset, fraction) {
    n.cells <- dim(dataset)[2]

    min.cells <- ceiling(fraction*n.cells)
    max.cells <- ceiling(fraction*n.cells)
    min.reads <- 2

    return(list(min.cells = min.cells, max.cells = max.cells,
                min.reads = min.reads))
}

#' Gene filter
#'
#' Filter genes that would not contribute to clustering, because they are either
#' expressed or not expressed in almost all cells.
#'
#' @param d Expression matrix with rows as genes and columns as cells
#' @param min.cells Minimum number of cells in which a given gene is expressed
#' (obtained from filter_params())
#' @param max.cells Maximum number of cells in which a given gene is expressed
#' (obtained from filter_params())
#' @param min.reads Minimum number of reads per gene per cell
#' (obtained from filter_params())
#' @return Filtered expression matrix in which only genes that are expressed in
#' more than \code{min.cells} with more than \code{min.reads} reads and also are
#' expressed in less than [total number of cells - \code{max.cells}].
gene_filter1 <- function(d, min.cells, max.cells, min.reads) {
    d <-
        d[rowSums(d > min.reads) >= min.cells &
              rowSums(d > 0) <= dim(d)[2] - max.cells, ]
    d <- unique(d)
    return(d)
}

#' Calculate a distance matrix
#'
#' Calculate a distance between column vectors of the input dataset using a
#' specified distance metrics.
#'
#' @param d Expression matrix with rows as genes and columns as cells
#' @param method Distance metrics: "spearman", "pearson", "euclidean", "maximum",
#' "manhattan", "canberra", "binary" or "minkowski"
#' @return A distance matrix
#'
#' @importFrom stats cor dist
#'
calculate_distance <- function(d, method) {
    return(if (method == "spearman") {
        as.matrix(1 - cor(d, method = "spearman"))
    } else if (method == "pearson") {
        as.matrix(1 - cor(d, method = "pearson"))
    } else {
        as.matrix(dist(t(d), method = method))
    })
}

#' Dimensionality reduction of a distance matrix
#'
#' Transform a distance matrix to a new basis
#'
#' @param dists A distance matrix
#' @param method Dimensionality reduction method: "pca", "spectral",
#' "spectral_reg" or "mds"
#' @return A transformed distance matrix
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

#' Laplacian
#'
#' Calculate Laplacian of an input distance matrix
#'
#' @param x Adjacency/distance matrix
#' @param tau Regularization term
#' @return Laplacian of the adjacency/distance matrix
norm_laplacian <- function(x, tau) {
    x <- x + tau * matrix(1, dim(x)[1], dim(x)[2])
    D <- diag(colSums(x))
    D1 <- D^(-0.5)
    D1[D1 == Inf] <- 0
    return(diag(dim(D)[1]) - D1 %*% x %*% D1)
}

#' Consensus clustering
#'
#' Caclulate consensus clustering by using Cluster-based similarity
#' partitioning algorithm (CSPA)
#'
#' @param clusts a list clustering labels (separated by a space)
#' @return Consensus matrix
#'
#' @importFrom stats dist
#'
consensus_clustering <- function(clusts) {
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

#' Support vector machines (SVM)
#'
#' Train an SVM classifier on teacher cells and then
#' classify study cells usin the classifier
#'
#' @param teach Expression matrix of teacher cells
#' @param study Expression matrix of study cells
#' @param kern Kernel to be used with SVM
#' @return Classification of study cells
#'
#' @importFrom e1071 svm
#' @importFrom stats predict
support_vector_machines <- function(teach, study, kern) {
    teach <- t(teach)
    labs <- factor(rownames(teach))
    rownames(teach) <- NULL
    model <- tryCatch(svm(teach, labs, kernel = kern),
                      error = function(cond) return(NA))
    pred <- predict(model, t(study))
    return(pred = pred)
}
