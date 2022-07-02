
#include <Rcpp.h>

// https://github.com/twolodzko/extraDistr/blob/master/src/beta-binomial-distribution.cpp

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;


// MACROS

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))


bool isInteger(double x, bool warn = true);


bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}


inline double logpmf_bbinom(double k, double n, double alpha,
                            double beta, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return k+n+alpha+beta;
#endif
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || !isInteger(n, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || k > n)
    return R_NegInf;
  // R::choose(n, k) * R::beta(k+alpha, n-k+beta) / R::beta(alpha, beta);
  return R::lchoose(n, k) + R::lbeta(k+alpha, n-k+beta) - R::lbeta(alpha, beta);
}


// https://github.com/twolodzko/extraDistr/blob/master/src/beta-binomial-distribution.cpp

// [[Rcpp::export]]
NumericVector cppdbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(),
                alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++){
    p[i] = logpmf_bbinom(GETV(x, i), GETV(size, i),
                         GETV(alpha, i), GETV(beta, i),
                         throw_warning);
  }

  if (!log_prob){
    p = Rcpp::exp(p);
  }
  
  if (throw_warning){
    Rcpp::warning("NaNs produced");
  }

  return p;
}
