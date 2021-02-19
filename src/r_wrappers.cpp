#include <RcppArmadillo.h>

#include "utility.hpp"
#include "rv.hpp"
#include "individual.hpp"
#include "cbga_r.hpp"

//[[Rcpp::export()]]
Rcpp::List cbga_proportional(
        Rcpp::Function & f, 
        const arma::colvec & lower, 
        const arma::colvec & upper,
        const int & max_generations, 
        const std::vector<int> & n_subpopulation, 
        const std::vector<int> n_chromosomes,
        const double & pi_mutation, 
        const double & pi_recombination, 
        const int & max_lifespan,
        const double & s_age, 
        const int & trace
) {
    //
    const int & n_genome = n_chromosomes.size();
    
    const int & n_total = std_vector_sum(
        n_chromosomes,
        n_genome
    );
    
    //
    RV random_variate;
    
    Fitness fitness(
            f,
            lower, 
            upper, 
            s_age
    );
    
    TotalPopulation donors(
            n_population[0], 
            n_genome, 
            n_chromosomes, 
            fitness, 
            random_variate
    );
    
    TotalPopulation breeders(
            n_population[1], 
            n_genome, 
            n_chromosomes, 
            fitness, 
            random_variate
    );
    
    TotalPopulation litter(
            n_population[1], 
            n_genome, 
            n_chromosomes, 
            fitness, 
            random_variate
    );
    
    //
    CBGAR CBGA(
            max_generations, 
            fitness, 
            pi_recombination, 
            pi_mutation,
            max_lifespan, 
            donors, 
            breeders, 
            litter,
            random_variate, 
            trace
    );
    
    CBGA.Run();
    
    //
    Rcpp::List RListReturn = CBGA.RListReturn();
    return RListReturn;
}

