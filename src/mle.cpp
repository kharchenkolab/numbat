#include <cmath>  // std::pow
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace roptim;

double l_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig, double phi = 1.0);

class fit_lnpois : public Functor {

    public:
    
    const std::vector<int> Y_obs;
    
    const std::vector<double> lambda_ref;

    const int d;

    // initialize with source and destination
    fit_lnpois(std::vector<int> Y_obs, const std::vector<double> lambda_ref, const int d): 
        Y_obs(Y_obs), lambda_ref(lambda_ref), d(d) {}

    double operator()(const arma::vec &x) override {
        return -l_lnpois_cpp(Y_obs, lambda_ref, d, x[0], x[1]);
    };
};

// [[Rcpp::export]]
arma::rowvec fit_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d) {

  fit_lnpois model(Y_obs, lambda_ref, d);

  Roptim<fit_lnpois> opt("L-BFGS-B");
  opt.control.trace = 0;
  opt.set_hessian(false);
  arma::vec lower = {-arma::datum::inf, 0.01};
  opt.set_lower(lower);

  arma::vec x = {0, 1};
  opt.minimize(model, x);

  return opt.par().t();
}

