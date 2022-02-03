
#include <Rcpp.h>

using namespace Rcpp;


Rcpp::List getjcpp(Rcpp::List x, int j) {

    //
    // Do in R if export:
    // if (x == NULL){
    //    return R_NilValue;
    //}
    //if (j <= 0) {  
    //    Rcpp::stop("j cannot be <= 0");
    //}

    int n = x.size();
    Rcpp::List copylist = clone(x);
    for (int i = 0; i < n; i++) {
        std::vector<int> elem = x[i];
        if (j > elem.size()) {  
            // j >= elem.size() + 1 due to 1-indexing in R
            Rcpp::stop("j larger than elements in the list");
        }
        copylist[i] = elem[j-1]; 
    }

    return copylist;
}


