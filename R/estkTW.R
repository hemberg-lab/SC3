#' Estimate the optimal k for k-means clustering
#' 
#' The procedure filters and log-transforms the expression matrix before 
#' finding the eigenvalues of the sample covariance matrix. It will then return 
#' the number of significant eigenvalues according to the Tracy-Widom test.
#' 
#' @param dataset processed input expression matrix.

estkTW <- function(dataset) {
  
  p <- ncol(dataset)
  n <- nrow(dataset)
  
  message("Computing eigenvalues...")
  
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
    if (evals[i] > bd) {
      k <- k + 1
    }
  }
  return(k)
}
