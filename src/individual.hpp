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
    double calculate_fitness(const arma::colvec & parameters);
    double calculate_fitness_scaled(const double & fitness, const double & age,
                                    const int & n_mutations);
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
    
    arma::colvec generate_random_chromosome(const int & N, RV & random_variate);
    
    // Constructors
    Individual();
    Individual(const int & _n_genome, const std::vector<int> & _n_chromosomes, 
               Fitness & _fitness, RV & _random_variate);
    
    double decode_chromosome(const arma::colvec & x, 
                             const double & l, const double & u, const int & N);
    void decode_individual(Fitness & fitness, 
                           const std::vector<int> & n_chromosomes, const int & n_genome);
};


class Population 
{
public:
    // Objects
    int n_population;
    int n_genome;
    std::vector<int> n_chromosomes;
    
    std::vector<Individual> p_active;
    std::vector<Individual> p_litter;
    double population_entropy;
    
    // Constructors
    Population(const int & _n_population, const int & _n_genome, 
               const std::vector<int> _n_chromosomes, 
               Fitness & _fitness, RV & _random_variate);
    
    void update_population_entropy();
};


#endif