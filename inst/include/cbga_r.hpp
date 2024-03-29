#ifndef cbga_r_hpp
#define cbga_r_hpp

#include <RcppArmadillo.h>

#include "rv.hpp"
#include "individual.hpp"

class CBGAR 
{ 
private:
    // Objects
    RV random_variate;
    Fitness fitness;
    Population population;
    
    double pi_recombination;
    double pi_mutation;
    
    int breeding_rotation_rate;
    arma::colvec breeding_rotation;
    
    int max_lifespan;
    int n_lifetime;
    int trace;
    
    // Parent selection
    int parent_selection(
            const arma::colvec & accumulated_fitness
    );
    
    // Recombination and mutation
    void recombine_mutate(
            Individual & C,
            const Individual & I1, 
            const Individual & I2
    );
    
    // Litter selection
    int litter_selection(
            const arma::colvec & accumulated_fitness
    );
    
public:
    // Objects
    Individual best_individual;
    
    // Constructors
    CBGAR(
        const int & _n_lifetime, 
        Fitness & _fitness,
        const double & _pi_recombination, 
        const double & _pi_mutation, 
        const int & _max_lifespan, 
        Population & _population,
        const int & _breeding_rotation_rate,
        const arma::colvec & _breeding_rotation,
        RV & _random_variate, 
        const int & _trace
    );
    
    // Functions
    void run();
    Rcpp::List RListReturn();
};

#endif