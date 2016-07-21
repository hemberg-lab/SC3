#include <sstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
std::vector<int> splits(const std::string myString) {

	std::stringstream iss(myString);

	int number;
	std::vector<int> myNumbers;
	while (iss >> number)
		myNumbers.push_back(number);
	return myNumbers;
}
// [[Rcpp::export]]
arma::mat calcPWD1(const arma::mat & x) {
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
// [[Rcpp::export]]
NumericMatrix calcPWD2(const NumericMatrix & x) {
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


// [[Rcpp::export]]
arma::mat consmx(const std::vector<std::string> myString, arma::mat res, int length) {

	std::vector<int> t;
	arma::mat s;
	arma::mat t1;

	int r;
	for (int i = 0; i < length; i += 1) {
		t = splits(myString[i]);
		t1 = conv_to<mat>::from(t);
		s = calcPWD1(t1);
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
