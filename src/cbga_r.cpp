#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "rv.hpp"
#include "individual.hpp"
#include "cbga_r.hpp"

// Constructor
CBGAR::CBGAR(const int & _n_lifetimes, Fitness & _fitness, 
             const double & _pi_recombination, const double & _pi_mutation,
             Population & _population, RV & _random_variate, const bool & _trace) :
    n_lifetimes(_n_lifetimes), fitness(_fitness), 
    pi_recombination(_pi_recombination), pi_mutation(_pi_mutation),
    population(_population), random_variate(_random_variate), trace(_trace) { }


void CBGAR::mutate() 
{
    for (int n = 0; n < population.n_population; n++) 
    {
        int n_mutations = 0;
        Individual & individual = population.p_litter[n];
        for (int i = 0; i < population.n_genome; i++) 
        {
            const int & n_chromosomes_i = population.n_chromosomes[i];
            for (int j = 0; j < n_chromosomes_i; j++) 
            {
                const double & u = random_variate.generate_uniform_real();
                if (u < pi_mutation) 
                {
                    
                    individual.chromosomes[i][j] = (1 - individual.chromosomes[i][j]);
                    n_mutations++;
                }
            }
        }
        
        individual.age = 0;
        individual.n_mutations = n_mutations;
        individual.decode_individual(fitness, population.n_chromosomes, population.n_genome);
        individual.fitness = fitness.calculate_fitness(individual.parameters);
        individual.fitness_scaled = fitness.calculate_fitness_scaled(individual.fitness, 
                                                                     individual.age, 
                                                                     individual.n_mutations);
        population.p_litter[n] = individual;
    }
}

// Public
void CBGAR::run()
{
    
}

Rcpp::List CBGAR::RListReturn() 
{
    return Rcpp::List::create(Rcpp::Named("X") = 0.0);
}
