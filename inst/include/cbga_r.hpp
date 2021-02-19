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
    
    TotalPopulation donors;
    TotalPopulation breeders;
    TotalPopulation litter;
    
    RV random_variate;
    
    double pi_recombination;
    double pi_mutation;
    
    int max_lifespan;
    int max_generations;
    int trace;
    
    // Selecting breeder for doner
    Individual breeder_selection(
            const arma::colvec & accumulated_fitness
    );
    
    // Recombination and mutation
    void recombine_mutate(
            Individual & C,
            const Individual & I1, 
            const Individual & I2
    );
    
    // Updating
    void update_donors();
    void update_breeders();
    void update_litter();
    
    
public:
    // Objects
    Individual super_individual;
    
    // Constructors
    CBGAR(
        const int & _max_generations, 
        Fitness & _fitness,
        const double & _pi_recombination, 
        const double & _pi_mutation, 
        const int & _max_lifespan, 
        TotalPopulation & _donors,
        TotalPopulation & _breeders,
        TotalPopulation & _litter,
        RV & _random_variate, 
        const int & _trace
    );
    
    // Functions
    void Run();
    Rcpp::List RListReturn();
};

#endif