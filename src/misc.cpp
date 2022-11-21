#include <Rcpp.h>
using namespace Rcpp;

// Based on https://github.com/cran/ape/blob/390386e67f9ff6cd8e6e523b7c43379a1551c565/src/plot_phylo.c
// [[Rcpp::export]]
NumericVector node_depth(int ntip, NumericVector e1, NumericVector e2,
        int nedge, NumericVector xx, int method)
/* method == 1: the node depths are proportional to the number of tips
   method == 2: the node depths are evenly spaced */
{

    int i;

    /* First set the coordinates for all tips */
    for (i = 0; i < ntip; i++) xx[i] = 1;

    /* Then compute recursively for the nodes; we assume `xx' has */
    /* been initialized with 0's which is true if it has been */
    /* created in R (the tree must be in pruningwise order) */
    if (method == 1) {
        for (i = 0; i < nedge; i++)
            xx[e1[i] - 1] = xx[e1[i] - 1] + xx[e2[i] - 1];
    } else { /* *method == 2 */
        for (i = 0; i < nedge; i++) {
            /* if a value > 0 has already been assigned to the ancestor
               node of this edge, check that the descendant node is not
               at the same level or more */
            if (xx[e1[i] - 1])
            if (xx[e1[i] - 1] >= xx[e2[i] - 1] + 1) continue;
            xx[e1[i] - 1] = xx[e2[i] - 1] + 1;
        }
    }
    return xx;
}

// gtools is orphaned
// [[Rcpp::export]]
int roman2int_internal(Rcpp::StringVector letters, int nchar) {

    std::vector<int> values;
    if (nchar < 1) {
        return NA_INTEGER;
    }

    for(int i=0; i<nchar; i++){
        if (letters[0][i]== 'I'){
            values.push_back(1);
        } else if(letters[0][i]== 'V'){
            values.push_back(5);
        } else if(letters[0][i]== 'X'){
            values.push_back(10);
        } else if(letters[0][i]== 'L'){
            values.push_back(50);
        } else if(letters[0][i]== 'C'){
            values.push_back(100);
        } else if(letters[0][i]== 'D'){
            values.push_back(500);
        } else if(letters[0][i]== 'M'){
            values.push_back(1000);
        } else {
            stop("Invalid roman numeral '%c'", letters[0][i]);
        }
    }

    int total = 0;
    if (nchar > 1) {
        for(int n=0; n<nchar-1; n++) {
            if(values[n] < values[n+1]){
                total-=values[n];
            } else {
                total+=values[n];
            }
        }
    }
    total += values[nchar-1];

    return total;
}