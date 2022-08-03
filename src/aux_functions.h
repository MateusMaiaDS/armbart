#include<iostream>
#include <vector>
#include <math.h>
#include<cmath>
#include<RcppArmadillo.h>
#include <algorithm>

using namespace std;

// Function to sample an integer from a sequence from
// 0:n
int sample_int(int n){
  return floor(R::runif(0,n));
}

// Sample an uniform value from a double vector
//[[Rcpp::export]]
double sample_double(arma::vec vec, int n_min_size){

  // Getting the range of the vector
  std::sort(vec.begin(),vec.end() );
  double value_index = floor(R::runif((n_min_size-1),(vec.size()-(n_min_size+1))));
  double value = vec(value_index);
  if(isnan(value)){
    cout << "ERROR" <<endl;
  }
  return value;

}


//[[Rcpp::export]]
arma::vec cpp_rgamma(int n, double shape, double scale){

  // Gamma samples
  arma::vec cpp_rgamma(n);

  for(int i=0;i<n;i++){
    cpp_rgamma(i) = R::rgamma(shape,scale);
  }
  return cpp_rgamma;

}
