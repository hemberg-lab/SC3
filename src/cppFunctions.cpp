#include <armadillo>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Compute Euclidean distance matrix by rows
//' 
//' Used in consmx function
//' 
//' @param x A numeric matrix.
// [[Rcpp::export]]
arma::mat ED1(const arma::mat & x) {
	unsigned int outrows = x.n_rows, i = 0, j = 0;
	double d;
	mat  out = zeros<mat>(outrows, outrows);

	for (i = 0; i < outrows - 1; i++) {
		arma::rowvec v1 = x.row(i);
		for (j = i + 1; j < outrows; j++) {
			d = sqrt(sum(pow(v1 - x.row(j), 2.0)));
			out(j, i) = d;
			out(i, j) = d;
		}
	}

	return out;
}

//' Compute Euclidean distance matrix by columns
//' 
//' Used in sc3-funcs.R distance matrix calculation
//' and within the consensus clustering.
//' 
//' @param x A numeric matrix.
// [[Rcpp::export]]
NumericMatrix ED2(const NumericMatrix & x) {
	unsigned int outcols = x.ncol(), i = 0, j = 0;
	double d;
	NumericMatrix out(outcols, outcols);

	for (j = 0; j < outcols - 1; j++) {
		NumericVector v1 = x.column(j);
		for (i = j + 1; i < outcols; i++) {
			d = sqrt(sum(pow(v1 - x.column(i), 2.0)));
			out(i, j) = d;
			out(j, i) = d;
		}
	}

	return out;
}

//' Consensus matrix computation
//' 
//' Computes consensus matrix given cluster labels
//' 
//' @param dat a matrix containing clustering solutions in columns
// [[Rcpp::export]]
arma::mat consmx(const arma::mat dat) {

	mat res = dat.n_cols * eye<mat>( dat.n_rows, dat.n_rows );

	int i, j, k;
	for (j = 0; j < dat.n_cols; j++) {
		for (i = 0; i < dat.n_rows; i++) {
			for (k = i + 1; k < dat.n_rows; k++) {
				if (dat(i, j) == dat(k, j)) {
				    res(i, k)++;
					res(k, i)++;
				}
			}
		}
	}
	res /= dat.n_cols;
	return res;
}

//' Graph Laplacian calculation
//' 
//' Calculate graph Laplacian of a symmetrix matrix
//' 
//' @param A symmetric matrix
//' @export
// [[Rcpp::export]]
arma::mat norm_laplacian(arma::mat A) {
    A = exp(-A/A.max());
    arma::mat D = diagmat(pow(sum(A), -0.5));
    arma::mat res = eye(A.n_cols, A.n_cols) - D * A * D;
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

