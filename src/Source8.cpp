#include <sstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

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