#include <Rcpp.h>
using namespace Rcpp;

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