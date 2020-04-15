#include <RcppArmadillo.h>

#include "utility.hpp"
#include "rv.hpp"
#include "individual.hpp"
#include "cbga_r.hpp"

//[[Rcpp::export()]]
Rcpp::List cbga_proportional(Rcpp::Function & f, const arma::colvec & lower, const arma::colvec & upper, 
                          const int & n_lifetimes, const int & n_population, const std::vector<int> n_chromosomes, 
                          const double & pi_mutation, const double & pi_recombination, 
                          const double & s_age, const double & s_mutation, 
                          const bool & trace) 
{
    //
    const int & n_genome = n_chromosomes.size();
    const int & n_total = std_vector_sum(n_chromosomes, n_genome);
    const double & e_mutation = pi_mutation * n_total;
    
    //
    RV random_variate;
    Fitness fitness(f, lower, upper, s_age, pi_mutation, e_mutation, s_mutation);
    Population population(n_population, n_genome, n_chromosomes, fitness, random_variate);
    
    //
    CBGAR CBGA(n_lifetimes, fitness, pi_recombination, pi_mutation,
               population, random_variate, trace);
    CBGA.run();
    
    //
    Rcpp::List RListReturn = CBGA.RListReturn();
    return RListReturn;
}