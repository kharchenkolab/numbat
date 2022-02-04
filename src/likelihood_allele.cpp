
#include <Rcpp.h>

using namespace Rcpp;


// expensive for loops in likelihood_allele() and forward_back_allele)()
// WIP


// [[Rcpp::export]]
Rcpp::NumericMatrix likelihood_allele_compute(Rcpp::List obj) {

    // calculate m, n
    Rcpp::NumericVector x = obj["x"]; // x <- obj$x
    Rcpp::List Pi = obj["Pi"];
    Rcpp::NumericMatrix Pi1 = Pi[0];
    int m = Pi1.nrow();  // m <- nrow(obj$Pi[[1]])
    int n = x.length();


    std::vector<double> objdelta =  Rcpp::as<std::vector<double>>(obj["delta"]);
    // std::valarray doesn't work with Rcpp?
    // std::for_each(objdelta.begin(), objdelta.end(), &std::log);
    for (int i = 0; i < objdelta.size(); i++) {
        objdelta[i] = std::log(objdelta[i]);
    }
    Rcpp::NumericVector logphi = Rcpp::wrap(objdelta); //logphi <- log(as.double(obj$delta))
    
    // logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
    const int nrow = n;
    const int ncol = m;
    Rcpp::NumericMatrix final(nrow, ncol); 

    return final;
}


