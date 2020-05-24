#ifndef cbga_r_hpp
#define cbga_r_hpp

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
    
    int n_lifetimes;
    bool trace;
    
    // Selection functions
    
    // Recombination
    Individual recombine(const int & n1, const int & n2);
    
    // Mutation
    void mutate();
    
    
public:
    // Objects
    Individual best_individual;
    
    // Constructors
    CBGAR(const int & _n_lifetimes, Fitness & _fitness, 
          const double & _pi_recombination, const double & _pi_mutation, 
          Population & _population, RV & _random_variate, const bool & _trace);
    
    // Functions
    void run();
    Rcpp::List RListReturn();
};

#endif