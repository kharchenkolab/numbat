#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppParallel;

// https://github.com/KlausVigo/phangorn/blob/master/src/phangorn_utils.cpp
// [[Rcpp::export]]
std::vector<std::vector<int>> allChildrenCPP(const arma::Mat<int> E) {

    arma::Col<int> parent = E.col(0);
    arma::Col<int> children = E.col(1);
    int m = max(parent);

    std::vector<std::vector<int>> out(m);

    for(int i = 0; i<parent.size(); i++) {
        out[parent(i)-1L].push_back(children(i));
    }

    return out;
}

// [[Rcpp::export]]
arma::mat CgetQ(arma::mat logQ, std::vector<std::vector<int>> children_dict, arma::Col<int> node_order){

    int n = node_order.n_rows;
    std::vector<int> children;

    for (int i = 0; i < n; ++i) {
        int node = node_order(i);
        children = children_dict[node-1];
        logQ.row(node-1) = logQ.row(children[0]-1) + logQ.row(children[1]-1);
    }

    return logQ;
}

// [[Rcpp::export]]
double score_tree_cpp(const arma::Mat<int> E, const arma::mat P) {

    int n = P.n_rows;
    int m = P.n_cols;

    arma::mat logQ(n * 2 - 1, m);

    arma::mat logP_0 = log(P);
    arma::mat logP_1 = log(1-P);

    logQ.rows(0, n-1) = logP_1 - logP_0;

    arma::Col<int> node_order(E.n_rows + 1);
    node_order.rows(0,E.n_rows-1) = E.col(1);
    node_order(E.n_rows) = n+1;
    arma::uvec ids = find(node_order > n);
    node_order = node_order.elem(ids);

    std::vector<std::vector<int>> children_dict = allChildrenCPP(E);

    logQ = CgetQ(logQ, children_dict, node_order);

    double l = 0;

    for (int i = 0; i < m; ++i) {
        l += max(logQ.col(i)) + sum(logP_0.col(i));
    }

    return l;
}

struct score_nni : public Worker {

    const std::vector<arma::Mat<int>> trees_vec;

    const arma::mat P;

    RVector<double> scores;

    // initialize with source and destination
    score_nni(const std::vector<arma::Mat<int>> trees_vec, const arma::mat P, NumericVector scores):
        trees_vec(trees_vec), P(P), scores(scores) {}

    // take the square root of the range of elements requested
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            scores[i] = score_tree_cpp(trees_vec[i], P);
        }
    }
};

// [[Rcpp::export]]
NumericVector score_nni_parallel(List trees, arma::mat P) {

    int n = trees.length();

    NumericVector scores(n);

    // Convert list to vector of double matrices
    std::vector<arma::Mat<int>> trees_vec;

    for (int i = 0; i < n; ++i) {
        List tree = trees[i];
        arma::Mat<int> E = tree["edge"];
        trees_vec.push_back(E);
    }

    score_nni score_nni(trees_vec, P, scores);

    parallelFor(0, n, score_nni);

    return scores;
}