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

/////////////////////////////////////// NNI ////////////////////////////////////////

// Modified from R-package `ape' by Emmanuel Paradis and Klaus Schliep

void bar_reorderRcpp(int node, int nTips, const arma::Col<int> & e1,
    const arma::Col<int> & e2, std::vector<int> & neworder, const arma::Col<int> & L,
    const arma::Col<int> & xi, const arma::Col<int> & xj, int & iii)
{
    int i = node - nTips - 1, j, k;

    for (j = xj[i] -1; j >= 0; j--)
        neworder[iii--] = L[xi[i] + j ] + 1;

    for (j = 0; j < xj[i]; j++) {
        k = e2[L[xi[i] + j ]];
        if (k > nTips)
            bar_reorderRcpp(k, nTips, e1, e2, neworder, L, xi, xj, iii);
    }
}

// [[Rcpp::export]]
arma::Mat<int> reorder_rows(arma::Mat<int> x, arma::Col<int> y) {

    // Create an output matrix
    arma::Mat<int> out = x;

    // Loop through each row and copy the data. 
    for (int i = 0; i < y.n_elem; ++i) {
        out.row(i) = x.row(y[i]-1);
    }

    return out;
}

// [[Rcpp::export]]
arma::Mat<int> reorderRcpp(arma::Mat<int> E) {

    int n = E.n_rows;
    int nTips = n/2 + 1;
    int root = nTips + 1;

    arma::Col<int> e1 = E.col(0);
    arma::Col<int> e2 = E.col(1);
    int m = max(e1), k, j;
    int nnode = m - nTips;
    
    arma::Col<int> L(n);
    std::vector<int> neworder(n);
    arma::Col<int> pos(nnode);
    arma::Col<int> xi(nnode);
    arma::Col<int> xj(nnode);
    for (int i = 0; i < n; i++) {
        xj[e1[i] - nTips - 1]++;
    }
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i-1] + xj[i - 1];
    }
    for (int i = 0; i < n; i++) {
        k = e1[i] - nTips - 1;
        j = pos[k]; /* the current 'column' position corresponding to k */
        L[xi[k] + j] = i;
        pos[k]++;
    }

    int iii = n - 1;

    bar_reorderRcpp(root, nTips, e1, e2, neworder, L, xi, xj, iii);

    E = reorder_rows(E, neworder);

    return E;
}

// Modified from R-package `phangorn' by Klaus Schliep
// [[Rcpp::export]]
std::vector<arma::Mat<int>> nnin_cpp(const arma::Mat<int> E, const int n) {

    arma::Mat<int> E1 = E;
    arma::Mat<int> E2 = E;
    arma::Col<int> parent = E.col(0);
    arma::Col<int> child = E.col(1);
    int k = min(parent) - 1;
    arma::uvec indvec = find(child > k);
    int ind = indvec[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    arma::uvec ind1_vec = find(parent == p1);
    ind1_vec = ind1_vec.elem(find(ind1_vec != ind));
    int ind1 = ind1_vec[0];
    arma::uvec ind2 = find(parent == p2);
    
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    E1(ind1, 1) = e2;
    E1(ind2[0], 1) = e1;
    E2(ind1, 1) = e3;
    E2(ind2[1], 1) = e1;

    std::vector<arma::Mat<int>> res(2);

    res[0] = reorderRcpp(E1);
    res[1] = reorderRcpp(E2);

    return res;
}

// Serial version
// [[Rcpp::export]]
List nni_cpp(const List tree) {
    
    arma::Mat<int> E = tree["edge"];
    int Nnode = tree["Nnode"];

    int n = E.n_rows/2 - 1;

    List res(2*n);

    for (int i = 0; i < n; i++) {
        std::vector<arma::Mat<int>> trees = nnin_cpp(E, i+1);
        arma::Mat<int> E1 = trees[0];
        arma::Mat<int> E2 = trees[1];
        
        List tree1 = List::create(Named("edge") = E1, Named("Nnode") = Nnode);
        List tree2 = List::create(Named("edge") = E2, Named("Nnode") = Nnode);

        res[2*i] = tree1;
        res[2*i+1] = tree2;
    }

    return res;

}

struct score_neighbours : public Worker {

    // original tree
    const arma::Mat<int> E;
    
    const arma::mat P;

    RVector<double> scores;

    // initialize with source and destination
    score_neighbours(const arma::Mat<int> E, const arma::mat P, NumericVector scores): 
        E(E), P(P), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Mat<int>> trees = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_cpp(trees[0], P);
            scores[2*i+1] = score_tree_cpp(trees[1], P);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel(const List tree, arma::mat P) {
    
    arma::Mat<int> E = tree["edge"];

    int n = E.n_rows/2 - 1;

    NumericVector scores(2*n);

    score_neighbours score_neighbours(E, P, scores);

    parallelFor(0, n, score_neighbours);

    return scores;

}
