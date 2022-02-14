
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

// [[Rcpp::export]]
double likelihood_allele_compute(Rcpp::List obj, Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, Rcpp::List logPi, int n, int m) {

    double LL = 0.0;

    for (int i = 0; i < n; i++) {

        if (i > 0) {
            Rcpp::NumericVector logphi_new; 
            for (int j = 0; j < m; j++) {
                Rcpp::NumericMatrix subset_logPi = logPi[i];    
                Rcpp::NumericVector logphi_logPi = logphi + subset_logPi(_, j);
                logphi_new.push_back(logSumExp(logphi_logPi));
            }
            logphi = logphi_new;
        }

        logphi = logphi + logprob(i, _);  // Note: logprob(i, _) is Rcpp::NumericVector

        double logSumPhi = logSumExp(logphi);

        logphi = logphi - logSumPhi;

        LL = LL + logSumPhi;
    }
    
    return LL;
}


// [[Rcpp::export]]
Rcpp::NumericVector forward_backward_compute(Rcpp::List obj, Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, Rcpp::List logPi, int n, int m) {

    const int nrow = n;
    const int ncol = m;
    Rcpp::NumericMatrix logalpha(nrow, ncol); // logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)

    double lscale = 0.0;
    double LL = 0.0;

    for (int t = 0; t < n; t++) {

        if (t > 0) {
            Rcpp::NumericVector logphi_new; 
            for (int j = 0; j < m; j++) {
                Rcpp::NumericMatrix subset_logPi = logPi[t];    
                Rcpp::NumericVector logphi_logPi = logphi + subset_logPi(_, j);
                logphi_new.push_back(logSumExp(logphi_logPi));
            }
            logphi = logphi_new;
        }

        logphi = logphi + logprob(t, _);  // Note: logprob(i, _) is Rcpp::NumericVector

        double logSumPhi = logSumExp(logphi);

        logphi = logphi - logSumPhi;

        lscale = lscale + logSumPhi;

        logalpha(t, _) = logphi + lscale;
    }

    LL = lscale;

    Rcpp::NumericMatrix logbeta(nrow, ncol); // logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
    Rcpp::NumericVector logphi_(m);  // logphi <- log(as.double(rep(1/m, m)))
    for (int i = 0; i < m; i++){
        double mval =  static_cast<double>(m);
        logphi_[i] = std::log(1/mval);
    }

    double lscale_ = std::log(m);

    for (int t = n-2; t>=0; t--) {

        // logphi = sapply(1:m, function(i) matrixStats::logSumExp(logphi + logprob[t+1,] + logPi[[t+1]][i,]))

        Rcpp::NumericVector logphinew;
        for (int i = 0; i < m; i++) {  
            Rcpp::NumericMatrix subset_logPi = logPi[t+1];   
            Rcpp::NumericVector logphi_logprob_logPi = logphi_ + logprob(t+1, _) + subset_logPi(i, _);
            logphinew.push_back(logSumExp(logphi_logprob_logPi));
        }
        logphi_ = logphinew;

        logbeta(t, _) = logphi_ + lscale_;
        
        double logSumPhi = logSumExp(logphi_);

        logphi_ = logphi_ - logSumPhi;

        lscale_ = lscale_ + logSumPhi;
        
    }

    // p_up = exp(logalpha + logbeta - LL)[,1]

    // matrix addition
    Rcpp::NumericMatrix expoutput(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            double sum = logalpha(i, j) + logbeta(i, j) - LL;
            expoutput(i, j) = std::exp(sum);  // Rcpp::exp() works with vectors
        }
    }

    // p_major/p_minor allele
    Rcpp::NumericVector p_major = expoutput(_, 1);

    return p_major;
}

