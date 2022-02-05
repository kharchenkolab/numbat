
#include <RcppArmadillo.h>

using namespace Rcpp;


#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
    // https://github.com/helske/seqHMM/blob/master/src/logSumExp.cpp
    unsigned int maxi = x.index_max();
    LDOUBLE maxv = x(maxi);
    if (!(maxv > -arma::datum::inf)) {
        return -arma::datum::inf;
    }
    LDOUBLE cumsum = 0.0;
    for (unsigned int i = 0; i < x.n_elem; i++) {
        if ((i != maxi) & (x(i) > -arma::datum::inf)) {
            cumsum += EXPL(x(i) - maxv);
        }
    }
  
    return maxv + log1p(cumsum);
}


// expensive for loops in likelihood_allele() and forward_back_allele)()
// WIP

// [[Rcpp::export]]
double likelihood_allele_compute(Rcpp::List obj, Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, Rcpp::List logPi, int n, int m) {

    // calculate m, n
    // Rcpp::NumericVector x = obj["x"]; // x <- obj$x
    //Rcpp::List Pi = obj["Pi"];
    //Rcpp::NumericMatrix Pi1 = Pi[0];
    //int m = Pi1.nrow();  // m <- nrow(obj$Pi[[1]])
    //int n = x.length();

    //std::vector<double> objdelta =  Rcpp::as<std::vector<double>>(obj["delta"]);

    //for (int i = 0; i < objdelta.size(); i++) {
    //    objdelta[i] = std::log(objdelta[i]);
    //}
    //Rcpp::NumericVector logphi = Rcpp::wrap(objdelta); //logphi <- log(as.double(obj$delta))
    
    //const int nrow = n;
    //const int ncol = m;
    // not used in function, Rcpp::NumericMatrix logalpha(nrow, ncol); // logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)

    double LL = 0.0;

    for (int i = 0; i < n; i++) {

        if (i > 1) {
            for (int j = 0; j < m; j++) {
                Rcpp::NumericVector subset_logPi = logPi[i];    
                Rcpp::NumericVector logphi_logPi = logphi + subset_logPi[j];
                logphi[j] = logSumExp(logphi_logPi);
            }
        }

        logphi = logphi + logprob(i, _);  // Note: logprob(i, _) is Rcpp::NumericVector

        double logSumPhi = logSumExp(logphi);

        logphi = logphi - logSumPhi;

        LL = LL + logSumPhi;
    }
    
    return LL;
}


