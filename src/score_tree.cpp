#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix CgetQ(NumericMatrix logQ, List children_dict, IntegerVector node_order){

    int n = node_order.size();
    IntegerVector children;

    for (int i = 0; i < n; ++i) {
        int node = node_order[i];
        children = children_dict[node-1];
        logQ(node-1,_) = logQ(children[0]-1,_) + logQ(children[1]-1,_);
    }

    return logQ;
}