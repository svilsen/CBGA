#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "rv.hpp"
#include "individual.hpp"
#include "cbga_r.hpp"
#include "utility.hpp"

// Constructor
CBGAR::CBGAR(
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
) :
    n_lifetime(_n_lifetime), fitness(_fitness), 
    pi_recombination(_pi_recombination), pi_mutation(_pi_mutation),
    max_lifespan(_max_lifespan), population(_population), 
    breeding_rotation_rate(_breeding_rotation_rate), breeding_rotation(_breeding_rotation),
    random_variate(_random_variate), trace(_trace) 
{ 
    best_individual = population.p_best[0];
    for (int k = 1; k < population.p_best.size(); k++) 
    {
        if (best_individual.fitness_scaled < population.p_best[k].fitness_scaled) 
        {
            best_individual = population.p_best[k];
        }
    }
}


int CBGAR::parent_selection(
        const arma::colvec & accumulated_fitness
) {
    const double & u = random_variate.generate_uniform_real();
    
    int n = 0;
    while (u > accumulated_fitness[n]) 
    {
        n++;
    }
    
    return n;
}

int CBGAR::litter_selection(
        const arma::colvec & accumulated_fitness
) {
    const double & u = random_variate.generate_uniform_real();
    
    int n = 0;
    while (u > accumulated_fitness[n]) 
    {
        n++;
    }
    
    return n;
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
    const int & K = population.n_population.size();
    for (int n = 0; n < n_lifetime; n++) 
    {
        for (int k = 0; k < K; k++) {
            // Update active accumulated proportional fitness
            arma::colvec accumulated_fitness_active_k = population.accumulated_proportional_fitness(
                population.p_active[k]
            );
            
            // Update litter population
            int M = population.n_population[k];
            for (int m = 0; m < M; m++) 
            {
                // Individual to override
                Individual C = population.p_litter[population.n_population_accumulate[k] + m];
                
                // Parent selection
                Individual I1 = population.p_best[k];
                
                int n = parent_selection(
                    accumulated_fitness_active_k
                );
                
                Individual I2 = population.p_active[k][n];
                
                // Recombination and mutation
                recombine_mutate(
                    C, 
                    I1, 
                    I2
                );
                
                // Create new decoded individual
                C.decode_individual(
                    fitness, 
                    population.n_chromosomes, 
                    population.n_genome
                );
                
                // Calculate new fitness
                C.fitness = fitness.calculate_fitness(
                    C.parameters
                );
                
                C.fitness_scaled = fitness.calculate_fitness_scaled(
                    C.fitness, 
                    C.age, 
                    C.n_mutations
                );
                
                // Override
                population.p_litter[population.n_population_accumulate[k] + m] = C;
            }
        }
        
        // Update active accumulated proportional fitness
        arma::colvec accumulated_fitness_litter = population.accumulated_proportional_fitness(
            population.p_litter
        );
        
        // Update active population, and best individual
        Individual active_best_individual;
        for (int k = 0; k < K; k++) {
            int M = population.n_population[k];
            for (int m = 0; m < M; m++) 
            {
                // Update active population
                Individual I = population.p_active[k][m];
                
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
                    int n = litter_selection(
                        accumulated_fitness_litter
                    );
                    
                    I = population.p_litter[n];
                }
                
                population.p_active[k][m] = I;
                
                
                // Update best individual in the active population
                if ((m == 0) || (I.fitness_scaled > active_best_individual.fitness_scaled)) {
                    active_best_individual = I;
                }
                
                // Update best individual in the entire lifetime
                if (I.fitness_scaled > population.p_best[k].fitness_scaled) 
                {
                    population.p_best[k] = I;
                }
            }
        }
        
        // Update population entropy
        population.update_population_entropy();
        
        // Update the stored best individual among the bulls, and update the 
        // ordering of p_best according to the breeding schedule.
        std::vector<Individual> p_best_new(K);
        for (int k = 0; k < K; k++) 
        {
            if (population.p_best[k].fitness_scaled > best_individual.fitness_scaled) 
            {
                best_individual = population.p_best[k];   
            }
            
            if (n % breeding_rotation_rate == 0) {
                p_best_new[k] = population.p_best[breeding_rotation[k]];
            }
        }
        
        if (n % breeding_rotation_rate == 0) {
            population.p_best = p_best_new;
        }
        
        if (trace > 0) 
        {
            double average_entropy = colvec_mean(population.population_entropy);
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
                            << "\t\tEntropy (average): " << average_entropy << "\n"
                            << "\t\tEntropy: " << population.population_entropy.t() << "\n";
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
                              Rcpp::Named("BreedingProtocol") = Rcpp::List::create(
                                  Rcpp::Named("Rate") = breeding_rotation_rate, 
                                  Rcpp::Named("Rotation") = breeding_rotation
                              ),
                              Rcpp::Named("ActiveEntropy") = population.population_entropy);
}
