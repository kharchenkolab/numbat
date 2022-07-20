#include <cmath>  // std::pow
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace roptim;
using namespace RcppParallel;
//using namespace Rcpp;


#include <Rcpp.h>

using namespace Rcpp;

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

struct fit_worker : public Worker {

    //const arma::Mat<int> count_mat;
    //const std::vector<double> lambda_ref;

    // input
    RMatrix<double> count_mat;
    RVector<double> lambda_ref;
    //std::vector<double> lambda_ref;

    // output matrix
    RMatrix<double> params;

    fit_worker(const Rcpp::NumericMatrix count_mat, Rcpp::NumericVector lambda_ref, Rcpp::NumericMatrix params): 
        count_mat(count_mat), lambda_ref(lambda_ref), params(params) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            //std::vector<int> counts = count_mat.column(i);
            //RcppParallel::RMatrix<int>::Column counts = count_mat.column(i);
            RcppParallel::RMatrix<double>::Column counts = count_mat.column(i);
            //int d = std::accumulate(counts.begin(), counts.end(), 0);
            int d = 2;
            //std::vector<int> counts_vec = static_cast<std::vector<int>>(counts);
            //std::vector<int> counts_vec = arma::conv_to< std::vector<int> >::from(counts);
            std::vector<int> counts_vec;
            for (int i = 0; i < counts.size(); i++) {
                counts_vec.push_back(counts[i]);
            }
            std::vector<double> lambdaRef;
             for (int i = 0; i < lambda_ref.size(); i++) {
                lambdaRef.push_back(lambda_ref[i]);
            }           
            arma::rowvec res = fit_lnpois_cpp(counts_vec, lambdaRef, d);
            params(i,0) = res(0);
            params(i,1) = res(1);
        }
    }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix fit_lnpois_parallel(Rcpp::NumericMatrix count_mat, Rcpp::NumericVector lambda_ref) {
    
    //int n = count_mat.n_cols;
    int n = 3;

    Rcpp::NumericMatrix params(n,2);

    fit_worker fit_worker(count_mat, lambda_ref, params);

    parallelFor(0, n, fit_worker);

    return params;

}


// class fit_phi : public Functor {

//     public:
    
//     const std::vector<int> Y_obs;
//     const std::vector<double> lambda_ref;
//     const int d;
//     const double mu;
//     const double sig;

//     // initialize with source and destination
//     fit_phi(std::vector<int> Y_obs, const std::vector<double> lambda_ref, const int d, const double mu, const double sig): 
//         Y_obs(Y_obs), lambda_ref(lambda_ref), d(d), mu(mu), sig(sig) {}

//     double operator()(const arma::vec &phi) override {
//         return -l_lnpois_cpp(Y_obs, lambda_ref, d, mu, sig, phi(0));
//     };
// };

// // [[Rcpp::export]]
// double fit_phi_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig) {

//   fit_phi model(Y_obs, lambda_ref, d, mu, sig);

//   Roptim<fit_phi> opt("L-BFGS-B");
//   opt.control.trace = 0;
//   opt.set_hessian(false);
//   arma::vec lower = {0.1};
//   arma::vec upper = {10};
//   opt.set_lower(lower);
//   opt.set_upper(upper);
//   arma::vec phi = {1};
//   opt.minimize(model, phi);
//   arma::vec phi_mle = opt.par();

//   return phi_mle(0);
// }