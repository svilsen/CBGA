#ifndef utility_hpp
#define utility_hpp

#include <RcppArmadillo.h>

int std_vector_sum(
        const std::vector<int> & x
);

std::vector<int> std_vector_accumulate(
        const std::vector<int> & x
);

double colvec_mean(
        const arma::colvec & x
);

#endif