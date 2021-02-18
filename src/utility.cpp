#include <Rcpp.h>
#include "utility.hpp"

int std_vector_sum(const std::vector<int> & x, const int & N) 
{
    int sum = 0.0;
    for (int n = 0; n < N; n++) 
    {
        sum += x[n];
    }
    
    return sum;
}
