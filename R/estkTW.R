#' Estimate the optimal k for k-means clustering
#' 
#' Based on random matrix theory
#' 
#' @param filename either an R matrix / data.frame object OR a
#' path to your input file containing an input expression matrix. The expression
#' matrix must contain both colnames (cell IDs) and rownames (gene IDs).
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
#

estkTW <- function(filename, 
                 gene.filter = TRUE,
                 gene.filter.fraction = 0.06,
                 gene.reads.rare = 2,
                 gene.reads.ubiq = 0,
                 log.scale = TRUE) {
  
  # get input data
  dataset <- get_data(filename)
  
  # remove duplicated genes
  dataset <- dataset[!duplicated(rownames(dataset)), ]
  
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
  
  p <- ncol(dataset)
  n <- nrow(dataset)
  cat("Read file with", n, "genes and", p, "cells.\n")
  
  cat("Computing eigenvalues...\n")
  
  # compute Tracy-Widom bound
  x <- scale(dataset)
  muTW <- (sqrt(n-1) + sqrt(p))^2
  sigmaTW <- (sqrt(n-1) + sqrt(p))*(1/sqrt(n-1) + 1/sqrt(p))^(1/3)
  sigmaHatNaive <- tmult(x) # x left-multiplied by its transpose
  bd <- 3.2730*sigmaTW + muTW # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
  
  # compute eigenvalues and return the amount which falls above the bound
  evals <- eigen(sigmaHatNaive, symmetric = TRUE, only.values = TRUE)$values
  k <- 0
  for (i in 1:length(evals)) {
    if (evals[i] > muTW) {
      k <- k + 1
    }
  }
  return(cat("Est. k =", k))
}