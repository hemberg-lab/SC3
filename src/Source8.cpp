#include <sstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <vector>
#include <string>
extern "C" void dsyevr_(char *jobz, char *range, char *uplo, int *n,
	double* a, int *lda, double *vl, double *vu,
	int *il, int *iu, double *abstol, int *m,
	double *w, double *z, int *ldz, int *isuppz,
	double *work, int *lwork, int *iwork,
	int *liwork, int *info);
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

// matrix multiplied by its transpose
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
	return(x.t()*x);
}


// [[Rcpp::export]]
int edec(Rcpp::NumericMatrix x, double bd = 0) {
  //double *A;
  int nrow = x.nrow();
  //A = new double[nrow*nrow];
  //int count = 0;
  //for (int i = 0; i < nrow; i += 1) {
    //for (int j = 0; j < nrow; j += 1) {
      //A[count] = x(i,j);
      //count += 1;
    //}
  //}
	char Jobz = 'N';
	char Range = 'V';
	char Uplo = 'U';
	int N = nrow;
	int LDA = N;
	double inf = std::numeric_limits<double>::infinity();
	int il;
	int iu;
	double abstol = 0.000001;
	int LDZ = N;
	int M;
	double W[N];
	double Z[N*N];
	int ISUPPZ[2*M];
	int LWORK = -1;
	int LIWORK = -1;
	double WWORK;
	int WIWORK;
	int INFO;

	dsyevr_(&Jobz, &Range, &Uplo, &N, &x(0,0), &LDA, &bd, &inf,
			&il, &iu, &abstol, &M, W, Z, &LDZ, ISUPPZ,
			&WWORK, &LWORK, &WIWORK, &LIWORK, &INFO);
	LWORK = WWORK;
	LIWORK = WIWORK;
	double WORK[LWORK];
	int IWORK[LIWORK];
	dsyevr_(&Jobz, &Range, &Uplo, &N, &x(0,0), &LDA, &bd, &inf,
         &il, &iu, &abstol, &M, W, Z, &LDZ, ISUPPZ,
         WORK, &LWORK, IWORK, &LIWORK, &INFO);
	return M;
}

// [[Rcpp::export]]
int ct(std::vector<double> evals, double bd, int length) {
  int k = 0;
  for (int i = 0; i < length; i += 1) {
    if (evals[i] > bd) {
      k += 1;
    }
  }
  return k;
}
