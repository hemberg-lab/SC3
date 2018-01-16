#include <RcppArmadillo.h>
#include <set>
#include <vector>
#include <algorithm>

using namespace arma;


//' Consensus matrix computation
//' 
//' Computes consensus matrix given cluster labels
//' 
//' @param dat a matrix containing clustering solutions in columns
//' @param K number of clusters
// [[Rcpp::export]]
arma::mat consmx(const arma::mat dat, int K) 
{
  using namespace std;
  typedef vector< set<int> > Membership;
  mat res = dat.n_cols * eye<mat>( dat.n_rows, dat.n_rows );
 
  vector <Membership> clusters;
 
  // Build index
  for (auto cm = 0; cm < dat.n_cols; cm++) 
  {
    Membership cluster(K , set<int>() );
		for (auto i = 0; i < dat.n_rows; i++) 
    {
      cluster[dat(i,cm)].insert(i);
		}
    // Push cluster lists to cluster membership container
    clusters.push_back(cluster);
  }
 
  // Build consensus matrix
  for (auto i1 = 0; i1 < clusters.size(); i1++)
  {
    // Cluster Result
    auto& cr1 = clusters[i1];
    for (auto i2 = i1 + 1; i2 < clusters.size(); i2++)
    {
      // Comparing clustering result
      auto& cr2 = clusters[i2];
      // Iterate through individual clusters in cr1
      for (auto& c1 : cr1)
      {
        for (auto& c2 : cr2)
        {
          vector <int> common_members;
          set_intersection(c1.begin(), c2.begin(), c1.end(), c2.end(), back_inserter(common_members));
          for (auto m1 = 0; m1 < common_members.size(); m1++)
          {
            for (auto m2 = m1 + 1; m2 < common_members.size(); m2++)
            {
              res(m1,m2)++;
              res(m2,m1)++;
            }
          }
        }
      }
    }
  }

  // Set diagonal back to one.. (Why? Nobody knows..)
  for ( auto i = 0; i < dat.n_rows; i++)
  {
    // is this legit??
    res(i,i) = 1;
    // check if armadillo support assignment operator
    // it should not matter that much, if not we can continue..
  }
  
  // Results should be now consistent with the original implementation, just faster
  
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
    arma::rowvec D_row = pow(sum(A), -0.5);
    A.each_row() %= D_row;
    arma::colvec D_col = conv_to< colvec >::from(D_row);
    A.each_col() %= D_col;
    arma::mat res = eye(A.n_cols, A.n_cols) - A;
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

