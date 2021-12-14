#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "rv.hpp"
#include "individual.hpp"
#include "cbga_r.hpp"

// Constructor
CBGAR::CBGAR(
    const int & _n_lifetime, 
    Fitness & _fitness,
    const double & _pi_recombination, 
    const double & _pi_mutation,
    const int & _max_lifespan, 
    Population & _population,
    RV & _random_variate, 
    const int & _trace
) :
    n_lifetime(_n_lifetime), fitness(_fitness), 
    pi_recombination(_pi_recombination), pi_mutation(_pi_mutation),
    max_lifespan(_max_lifespan), population(_population), 
    random_variate(_random_variate), trace(_trace) 
{
    best_individual = population.p_active[0];
    const int & M = population.n_population;
    for (int m = 1; m < M; m++)
    {
        Individual I = population.p_active[m];
        if (I.fitness > best_individual.fitness) 
        {
            best_individual = I;
        }
    }
}


Individual CBGAR::parent_selection(
        const arma::colvec & accumulated_fitness
) {
    const double & u = random_variate.generate_uniform_real();
    
    int n = 0;
    while (u > accumulated_fitness[n]) 
    {
        n++;
    }
    
    return population.p_active[n];
}

Individual CBGAR::litter_selection(
        const arma::colvec & accumulated_fitness
) {
    const double & u = random_variate.generate_uniform_real();
    
    int n = 0;
    while (u > accumulated_fitness[n]) 
    {
        n++;
    }
    
    return population.p_litter[n];
}

void CBGAR::recombine_mutate(
        Individual & C, 
        const Individual & I1, 
        const Individual & I2
) {
    int k = 0;
    std::vector<Individual> IA{I1, I2};
    if (I1.fitness_scaled < I2.fitness_scaled) 
    {
        k = 1;
    }
    
    int n_mutations = 0;
    for (int i = 0; i < population.n_genome; i++) 
    {
        const int & n_chromosomes_i = population.n_chromosomes[i];
        for (int j = 0; j < n_chromosomes_i; j++) 
        {
            // Recombination
            double u = random_variate.generate_uniform_real();
            if (u < pi_recombination) 
            {
                k = (k + 1) % 2;
            }
            
            C.chromosomes[i][j] = IA[k].chromosomes[i][j];
            
            // Mutation
            u = random_variate.generate_uniform_real();
            if (u < pi_mutation) 
            {
                C.chromosomes[i][j] = (1 - C.chromosomes[i][j]);
                n_mutations++;
            }
        }
    }
    
    C.age = 0;
    C.n_mutations = n_mutations;
}


// Public
void CBGAR::run()
{
    const int & M = population.n_population;
    for (int n = 0; n < n_lifetime; n++) 
    {
        // Update active accumulated proportional fitness
        arma::colvec accumulated_fitness_active = population.accumulated_proportional_fitness(
            population.p_active
        );
        
        // Update litter population
        for (int m = 0; m < M; m++) 
        {
            //
            Individual C = population.p_litter[m];
            
            //
            Individual I1 = best_individual;
            
            Individual I2 = parent_selection(
                accumulated_fitness_active
            );
            
            //
            recombine_mutate(
                C, 
                I1, 
                I2
            );

            //
            C.decode_individual(
                fitness, 
                population.n_chromosomes, 
                population.n_genome
            );
            
            C.fitness = fitness.calculate_fitness(
                C.parameters
            );
            
            C.fitness_scaled = fitness.calculate_fitness_scaled(
                C.fitness, 
                C.age, 
                C.n_mutations
            );
            
            population.p_litter[m] = C;
        }
        
        // Update active accumulated proportional fitness
        arma::colvec accumulated_fitness_litter = population.accumulated_proportional_fitness(
            population.p_litter
        );
        
        // Update active population, and best individual
        Individual active_best_individual;
        for (int m = 0; m < M; m++) 
        {
            // Update active population
            Individual I = population.p_active[m];
            
            double new_age = (max_lifespan * I.age + 1.0) / max_lifespan;
            double pi_new_age = 4.0 * (-(new_age * (1.0 - new_age)) + 0.25);
            if (new_age <= 0.5) 
            {
                pi_new_age = pi_new_age / 4;
            }
            
            double u = random_variate.generate_uniform_real();
            if (u > pi_new_age) 
            {
                I.age = new_age;
                I.fitness_scaled = fitness.calculate_fitness_scaled(
                    I.fitness, 
                    I.age, 
                    I.n_mutations
                );
            }
            else 
            {
                I = litter_selection(
                    accumulated_fitness_litter
                );
            }
            
            population.p_active[m] = I;
            
            
            // Update best individual in the active population
            if ((m == 0) || (I.fitness_scaled > active_best_individual.fitness_scaled)) {
                active_best_individual = I;
            }
            
            // Update best individual in the entire lifetime
            if (I.fitness_scaled > best_individual.fitness_scaled) 
            {
                best_individual = I;
            }
        }
        
        population.update_population_entropy();
        
        if (trace > 0) 
        {
            if ((n == 0) || ((n + 1) % trace) == 0) 
            {
                Rcpp::Rcout << "Iteration: " << n + 1 << " / " << n_lifetime << "\n"
                            << "\tBest individual (total):\n"
                            << "\t\tParameters: " << best_individual.parameters.t()
                            << "\t\tFitness: " << best_individual.fitness << "\n"
                            << "\t\tScaled fitness: " << best_individual.fitness_scaled << "\n"
                            << "\tBest individual (active):\n"
                            << "\t\tParameters: " << active_best_individual.parameters.t()
                            << "\t\tFitness: " << active_best_individual.fitness << "\n"
                            << "\t\tScaled fitness: " << active_best_individual.fitness_scaled << "\n"
                            << "\tActive population:\n"
                            << "\t\tEntropy: " << population.population_entropy << "\n";
            }
        }
    }
}

Rcpp::List CBGAR::RListReturn() 
{
    return Rcpp::List::create(Rcpp::Named("chromosomes") = best_individual.chromosomes, 
                              Rcpp::Named("parameters") = best_individual.parameters, 
                              Rcpp::Named("mutations") = best_individual.n_mutations, 
                              Rcpp::Named("age") = best_individual.age, 
                              Rcpp::Named("f") = best_individual.fitness, 
                              Rcpp::Named("fs") = best_individual.fitness_scaled, 
                              Rcpp::Named("ActiveEntropy") = population.population_entropy);
}
