#include <sstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

// Used for Spectral transformation
// More efficient way to pre and post multiply by diagonal matrix D
// [[Rcpp::export]]
arma::mat mult(arma::mat D, arma::mat x, int dim) {
	arma::mat res(dim, dim);
	for (int i = 0; i < dim; i += 1) {
		for (int j = 0; j < dim; j += 1) {
			res(i, j) = x(i, j)*D(i, i)*D(j, j);
		}
	}
	return(res);
}

// matrix left-multiplied by its transpose
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
	return(x.t()*x);
}



