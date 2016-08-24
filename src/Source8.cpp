#include <sstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

//' More efficient way to pre and post multiply by diagonal matrix D
//' 
//' If D is a diagonal matrix and A another matrix with the same
//' dimensions, this procedure is a more efficient way to compute
//' D %*% A %*% D.
//' 
//' @param D Diagonal numeric matrix.
//' @param x Numeric matrix with the same dimensions as D.
//' @param dim Number of rows of D.
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

//' Matrix left-multiplied by its transpose
//' 
//' Given matrix A, the procedure returns A'A.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
	return(x.t()*x);
}



