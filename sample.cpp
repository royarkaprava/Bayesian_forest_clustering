#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// arma::uword samp_temp(const arma::vec &pvec) {
//   
//   // generate a vector of indices, so that 0 represents the largest 
//   
//   // element, 1 the second largest, and so on
//   
//   arma::uvec indx = arma::sort_index(pvec, "descend");
//   
//   // generate the q-vector, the vector of cumulative sums
//   
//   arma::vec qvec = arma::cumsum(arma::sort(pvec, "descend"));
//   
//   // draw randomly from (0,1)
//   
//   double u = arma::randu();
//   
//   // find interval into which u falls
//   
//   for (arma::uword k = 0; k < pvec.n_elem; ++k) {
//     
//     if (u <= qvec(k))
//       return indx(k);
//     
//   }
// }

// [[Rcpp::export]]
void draw(const int n, const int k, mat &res) {
  
  res.row(0) = conv_to<rowvec>::from(resi);
  

}