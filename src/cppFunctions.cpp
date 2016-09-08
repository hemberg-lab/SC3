#include <sstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;
using namespace Rcpp;

//' Split integral string into integral vector 
//' 
//' Similar to "strsplit" in R but only for integral strings and vectors. 
//' 
//' To remove a note from check change this function using an answer here:
//' http://stackoverflow.com/questions/5607589/right-way-to-split-an-stdstring-into-a-vectorstring
//' 
//' @param myString An integral string.
// [[Rcpp::export]]

std::vector<int> splits(const std::string myString) {

	std::stringstream iss(myString);

	int number;
	std::vector<int> myNumbers;
	while (iss >> number)
		myNumbers.push_back(number);
	return myNumbers;
}

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

//' Consensus binary matrix computation
//' 
//' Computes binary matrix given cluster labels
//' 
//' @param res A square zero matrix, nrows = length. Will become
//' binary matrix.
//' @param myString Integral string corresponding to cluster labels.
//' @param length Integer, number of cells.
// [[Rcpp::export]]
arma::mat consmx(const std::vector<std::string> myString, arma::mat res, int length) {

	std::vector<int> t;
	arma::mat s;
	arma::mat t1;

	int r;
	for (int i = 0; i < length; i += 1) {
		t = splits(myString[i]);
		t1 = conv_to<mat>::from(t);
		s = ED1(t1);
		r = s.n_rows;
		for (int j = 0; j < r; j += 1) {
			for (int k = j; k < r; k += 1) { //we make use of the symmetry of the dist mx
				if (s(j, k) == 0) {
					s(j, k) = -1;
					s(k, j) = -1;
					//if (j != k) {
					///s(k, j) = -1;
					//}
				}
				if (s(j, k) != -1) {
					s(j, k) = 0;
					s(k, j) = 0;
					//if (j != k) {
					//s(k, j) = 0;
					//}
				}
			}
		}
		res += s;
	}
	return res;
}

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

