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
  mat res = zeros(dat.n_rows, dat.n_rows);
 
  vector <Membership> clusters;
 
  // Build index
  for (size_t cm = 0; cm < dat.n_cols; cm++) 
  {
    Membership cluster(K , set<int>() );
    for (size_t i = 0; i < dat.n_rows; i++) 
    {
      cluster[dat(i,cm) - 1].insert(i);
    }
    
    // Push cluster lists to cluster membership container
    clusters.push_back(cluster);
  }
 
  // Build consensus matrix
  for (size_t i1 = 0; i1 < clusters.size(); i1++)
  {
    
    // Cluster Result
    Membership& cr1 = clusters[i1];
    for (size_t i2 = i1 + 1; i2 < clusters.size(); i2++)
    {
      // Comparing clustering result
      Membership& cr2 = clusters[i2];
      // Iterate through individual clusters in cr1
      for (Membership::const_iterator c1 = cr1.begin();  c1 != cr1.end(); ++c1)
      {
        
        for (Membership::const_iterator c2 = cr2.begin();  c2 != cr2.end(); ++c2)
        {
          vector <int> common_members;
          set_intersection(c1->begin(), c1->end(), c2->begin(), c2->end(), back_inserter(common_members));
          for (size_t m1 = 0; m1 < common_members.size(); m1++)
          {
            for (size_t m2 = m1 + 1; m2 < common_members.size(); m2++)
            {
              res(common_members[m1], common_members[m2])++;
              res(common_members[m2], common_members[m1])++;
            }
          }
        }
      }
    }
  }

  // Set diagonal back to one.. (Why? Nobody knows..)
  for ( size_t i = 0; i < dat.n_rows; i++)
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

//' Matrix left-multiplied by its transpose
//' 
//' Given matrix A, the procedure returns A'A.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
    return(x.t()*x);
}

//' Converts the distance matrix to adjacency matrix
//' 
//' Given matrix A, the procedure returns a transformed matrix A'.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat distance_to_adjacency_mat(arma::mat A) {
    return exp(A*A.max());
}


