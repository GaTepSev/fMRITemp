#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector data_generate_cpp(List theta, List covariates, int N, int S, int T, int setting, NumericVector epsilon){
  double sigma2_nu = theta["sigma2_nu"];
  double sigma2_xi = theta["sigma2_xi"];
  NumericVector beta = theta["beta"];
  NumericVector x = covariates["x"];

  NumericVector nu = rnorm(N, 0.0, sigma2_nu);
  NumericVector xi = rnorm(N * S * T, 0.0, sigma2_xi);
  NumericVector y(N * S * T);

  int i, s, t, c = 0;
  for(i = 0; i < N; i++)
  {
    for(s = 0; s < S; s++)
    {
      for(t = 0; t < T; t++)
      {
        y[i + s * N + t * N * S] = x[i + N * t] * beta[s] + nu[i] + epsilon[s + t * S] + xi[i + s * N + t * N * S];
      }
    }
  }

  return y;
}
