
// https://github.com/KlausVigo/phangorn/blob/master/src/phangorn_utils.cpp

#include <Rcpp.h>

using namespace Rcpp;


// shorter and easier to understand replacement of C function
// import: edge matrix
// export: list of children


// [[Rcpp::export]]
List allChildrenCPP(const IntegerMatrix orig) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent);
    // create list for results
    std::vector< std::vector<int> > out(m) ;
    for(int i = 0; i<parent.size(); i++){
        out[parent[i]-1L].push_back(children[i]);
    }
    return wrap(out);
}
