#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "rv.hpp"
#include "individual.hpp"
#include "cbga_r.hpp"

// Constructor
CBGAR::CBGAR(const Fitness & _fitness, 
             const double & _pi_recombination, const double & _pi_mutation, 
             Population & _population, RV & _random_variate, const bool & _trace) :
    fitness(_fitness), 
    pi_recombination(_pi_recombination), pi_mutation(_pi_mutation),
    population(_population), random_variate(_random_variate), trace(_trace) { }


// Public


