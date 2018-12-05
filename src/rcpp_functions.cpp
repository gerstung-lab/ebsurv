#include <Rcpp.h>
using namespace Rcpp;

//' Multiply a number by two
//' 
//' @param time_vector,integrand_1,integrand_2  Numeric vectors.
//' @export
// [[Rcpp::export]]
NumericVector convolute_semiMarkov(NumericVector time_vector, NumericVector integrand_1, NumericVector integrand_2 ) {
  int xLen = time_vector.size();
  NumericVector overallSurvival(xLen);
  for(int i = 0; i < xLen; ++i) overallSurvival[i] = 0;
  for(int i = 1; i < xLen; ++i){
    for(int j = 1; j < i; ++j){
      overallSurvival[i] += integrand_2[i-j]*integrand_1[j];
    }
  }
  return overallSurvival;
}
