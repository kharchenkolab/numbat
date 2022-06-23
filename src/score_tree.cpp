// #include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppParallel;

// https://github.com/KlausVigo/phangorn/blob/master/src/phangorn_utils.cpp
// [[Rcpp::export]]
List allChildrenCPP(const IntegerMatrix orig) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent);
    // create list for results
    std::vector< std::vector<int> > out(m);
    for(int i = 0; i<parent.size(); i++){
        out[parent[i]-1L].push_back(children[i]);
    }
    return wrap(out);
}

// [[Rcpp::export]]
arma::mat CgetQ(arma::mat logQ, List children_dict, IntegerVector node_order){

    int n = node_order.size();
    IntegerVector children;

    for (int i = 0; i < n; ++i) {
        int node = node_order[i];
        children = children_dict[node-1];
        logQ.row(node-1) = logQ.row(children[0]-1) + logQ.row(children[1]-1);
    }

    return logQ;
}

// [[Rcpp::export]]
double score_tree_cpp(List tree, arma::mat P) {

    int k = tree["Nnode"];
    int n = P.n_rows;
    int m = P.n_cols;

    arma::mat logQ(k * 2 + 1, m);

    arma::mat logP_0 = log(P);
    arma::mat logP_1 = log(1-P);

    logQ.rows(0, n-1) = logP_1 - logP_0;

    IntegerMatrix E = tree["edge"];

    IntegerVector node_order = E(_,1);
    node_order.push_back(n+1);
    node_order = node_order[node_order > n];

    // List children_dict = allChildrenCPP(E);

    // logQ = CgetQ(logQ, children_dict, node_order);

    double l_tree = 0;

    // for (int i = 0; i < m; ++i) {
    //     l_tree += max(logQ.col(i)) + sum(logP_0.col(i));
    // }

    return l_tree;
}

struct score_nni : public Worker {

    const List trees;

    const arma::mat P;

    NumericVector scores;

    // initialize with source and destination
    score_nni(const List trees, arma::mat P, NumericVector scores):
        trees(trees), P(P), scores(scores) {}

    // take the square root of the range of elements requested
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            List tree = trees[i];
            // scores[i] = tree["Nnode"];
            // scores[i] = score_tree_cpp(tree, P);
            scores[i] = 0;
        }
    }
};

// [[Rcpp::export]]
NumericVector score_nni_parallel(List trees, arma::mat P) {

    int n = trees.length();

    NumericVector scores(n);

    score_nni score_nni(trees, P, scores);

    parallelFor(0, n, score_nni);

    // for (int i = 0; i < n; ++i) {
    //     List tree = trees[i];
    //     scores[i] = score_tree_cpp(tree, P);
    // }

    return scores;
}