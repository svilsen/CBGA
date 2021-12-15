#include <RcppArmadillo.h>
#include "utility.hpp"

int std_vector_sum(
        const std::vector<int> & x
) {
    const int & N = x.size(); 
    int sum = 0.0;
    for (int n = 0; n < N; n++) 
    {
        sum += x[n];
    }
    
    return sum;
}

std::vector<int> std_vector_accumulate(
        const std::vector<int> & x
) {
    const int & N = x.size(); 
    std::vector<int> xa(N + 1); 
    xa[0] = 0;
    for (int n = 0; n < N; n++) 
    {
        xa[n + 1] = xa[n] + x[n];
    }
    
    return xa;
}

double colvec_mean(
        const arma::colvec & x
) {
    const int & N = x.size(); 
    
    double x_sum = 0.0;
    for (int n = 0; n < N; n++) 
    {
        x_sum += x[n];
    }
    
    return x_sum / N;
}