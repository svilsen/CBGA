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
    
    double pi_mutation;
    double e_mutation;
    double s_mutation;
    
    // Functions
    double age_scaling(const double & age);
    double mutation_scaling(const int & n_mutations);
    
public:
    // Objects
    arma::colvec lower;
    arma::colvec upper;
    
    // Constructor
    Fitness(Rcpp::Function & _f, const arma::colvec & _lower, const arma::colvec & _upper, 
            const double & _s_age, const double & _pi_mutation, const double & _e_mutation, 
            const double & _s_mutation);
    
    // Functions
    double decode_chromosome(const arma::colvec & x, const double & l, const double & u, const int & N);
    void decode_individual(Individual & individual, const std::vector<int> & n_chromosomes);
    
    double calculate_fitness(const arma::colvec & paramters, const double & age, const int & n_mutations);
};

class Population 
{
private:
    // Objects
    int n_population;
    int n_genome;
    std::vector<int> n_chromosomes;
    
public:
    // Objects
    std::vector<Individual> p_active;
    std::vector<Individual> p_litter;
    double inbreeding_coefficient;
    
    // Constructors
    population(const int & _n_population, const int & _n_genome, const std::vector<int> _n_chromosome, 
               const Fitness & _fitness, const RV & _random_variate);
    
    void update_inbreeding_coefficient();
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
    
    arma::colvec generate_random_chromosome(const int & N);
    
    // Constructors
    Individual();
    Individual(const int & _n_genome, const std::vector<int> & _n_chromosomes, 
               const Fitness & _fitness, const RV & _random_variate);
};

#endif