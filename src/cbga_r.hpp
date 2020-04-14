#ifndef cbga_r_hpp
#define cbga_r_hpp

#include <Rcpp.h>
#include <RcppArmadillo.h>

#include "rv.hpp"
#include "individual.hpp"

class CBGAR 
{ 
private:
    // Objects
    Fitness fitness;
    
    Population population;
    RV random_variate;
    
    double pi_recombination;
    double pi_mutation;
    
    bool trace;
    
    // Utility functions
    void update_fitness(const int & n);
    
    // Selection functions
    
    // Mutation
    
    // Recombination
    
    
public:
    // Objects
    individual best_individual;
    
    // Constructors
    CBGAR(const Fitness & _fitness, 
          const double & _pi_recombination, const double & _pi_mutation, 
          Population & _population, const RV & _random_variate, const bool & _trace;
    
    // Functions
    void run();
};

#endif