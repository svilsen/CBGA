#ifndef individual_hpp
#define individual_hpp

#include <RcppArmadillo.h>
#include "rv.hpp"

class Fitness 
{
private:
    // Objects
    Rcpp::Function f;
    
    double s_age;
    
    // Functions
    double age_scaling(
            const double & age
    );
    
public:
    // Objects
    arma::colvec lower;
    arma::colvec upper;
    
    // Constructor
    Fitness(
        Rcpp::Function & _f, 
        const arma::colvec & _lower, 
        const arma::colvec & _upper, 
        const double & _s_age
    );
    
    // Functions
    double calculate_fitness(
            const arma::colvec & parameters
    );
    
    double calculate_fitness_scaled(
            const double & fitness, 
            const double & age
    );
};


class Individual 
{ 
public:
    // Objects
    std::vector<arma::colvec> chromosomes;
    arma::colvec parameters;
    
    int n_mutations;
    double age;
    
    double fitness;
    double fitness_scaled;
    
    // Functions
    arma::colvec generate_random_chromosome(
            const int & N, RV & random_variate
    );
    
    // Constructors
    Individual();
    Individual(
        const int & _n_genome, 
        const std::vector<int> & _n_chromosomes,
        Fitness & _fitness, 
        RV & _random_variate
    );
    
    // Functions
    double decode_chromosome(
            const arma::colvec & x,
            const double & l, 
            const double & u, 
            const int & N
    );
    void decode_individual(
            Fitness & fitness,
            const std::vector<int> & n_chromosomes, 
            const int & n_genome
    );
};


class SubPopulation 
{
public:
    // Objects
    int n_subpopulation;
    int n_genome;
    std::vector<int> n_chromosomes;
    
    std::vector<Individual> individuals;
    double population_entropy;
    
    // Constructors
    SubPopulation();
    SubPopulation(
        const int & _n_subpopulation, 
        const int & _n_genome, 
        const std::vector<int> _n_chromosomes, 
        Fitness & _fitness, RV & _random_variate
    );
    
    // Functions
    void update_population_entropy();
    arma::colvec accumulated_proportional_fitness(
            const std::vector<Individual> & p
    );
};

class TotalPopulation 
{
public:
    // Objects
    int n_population;
    
    std::vector<SubPopulation> populations;
    double average_population_entropy;
    
    // Constructors
    TotalPopulation(
        const int & _n_population, 
        const int & _n_genome, 
        const std::vector<int> _n_chromosomes, 
        Fitness & _fitness, RV & _random_variate
    );
    
    void update_population_entropy();
};


#endif