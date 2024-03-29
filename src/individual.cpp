#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "rv.hpp"
#include "individual.hpp"
#include "utility.hpp"

//// Fitness
// Constructor
Fitness::Fitness(
    Rcpp::Function & _f, 
    const arma::colvec & _lower, 
    const arma::colvec & _upper, 
    const double & _s_age, 
    const double & _pi_mutation, 
    const double & _e_mutation,
    const double & _s_mutation
) : 
    f(_f), 
    lower(_lower), 
    upper(_upper), 
    s_age(_s_age), 
    pi_mutation(_pi_mutation),
    e_mutation(_e_mutation),
    s_mutation(_s_mutation) 
{ }

// Functions
double Fitness::age_scaling(
        const double & age
) {
    double age_scaled = s_age * age * (1.0 - age);
    if ((age_scaled > 1.0) || (s_age == 0)) 
    {
        age_scaled = 1.0;
    }
    
    return age_scaled;
}

double Fitness::mutation_scaling(
        const int & n_mutations
) {
    double v = std::pow(pi_mutation * (1.0 - pi_mutation), 0.5);
    double r_standard = (n_mutations - e_mutation) / v;
    
    return std::exp(-s_mutation * std::abs(r_standard));
}

double Fitness::calculate_fitness(
        const arma::colvec & parameters
) {
    double fitness = Rcpp::as<double>(f(parameters));
    return fitness;
}

double Fitness::calculate_fitness_scaled(
        const double & fitness, 
        const double & age, 
        const int & n_mutations
) {
    double age_s = age_scaling(age);
    double mutation_s = mutation_scaling(n_mutations);
    
    return age_s * mutation_s * fitness;
}


//// Population 
// Constructor
Population::Population(
    const std::vector<int> & _n_population, 
    const int & _n_genome,
    const std::vector<int> _n_chromosomes,
    Fitness & _fitness, 
    RV & _random_variate
) :
    n_population(_n_population), 
    n_genome(_n_genome),
    n_chromosomes(_n_chromosomes)
{
    //
    n_population_accumulate = std_vector_accumulate(n_population); 
    
    int K = n_population.size();
    int n_population_sum = n_population_accumulate[K];
    
    p_best = std::vector<Individual>(K);
    p_active = std::vector<std::vector<Individual>>(K);
    p_litter = std::vector<Individual>(n_population_sum);
    
    for (int k = 0; k < K; k++) 
    {
        int N = n_population[k];
        Individual p_best_k;
        std::vector<Individual> p_active_k = std::vector<Individual>(N);
        for (int n = 0; n < N; n++) 
        {
            Individual individual(
                    n_genome, 
                    n_chromosomes, 
                    _fitness, 
                    _random_variate
            );
            
            p_active_k[n] = individual;
            p_litter[n_population_accumulate[k] + n] = individual;
            
            if ((n == 0) || (individual.fitness_scaled > p_best_k.fitness_scaled)) {
                p_best_k = individual;
            }
        }
        
        p_active[k] = p_active_k;
        p_best[k] = p_best_k;
    }
    
    update_population_entropy();
}

void Population::update_population_entropy()
{
    population_entropy = arma::colvec(n_population.size());
    for (int k = 0; k < n_population.size(); k++) 
    {
        double population_entropy_k = 0.0;
        for (int i = 0; i < n_genome; i++) 
        {
            const int & n_chromosomes_i = n_chromosomes[i];
            for (int j = 0; j < n_chromosomes_i; j++)
            {
                double n_kij = 1.0;
                for (int n = 0; n < n_population[k]; n++) 
                {
                    n_kij += p_active[k][n].chromosomes[i][j];
                }
                
                const double p_kij = n_kij / (n_population[k] + 1.0);
                const double e_kij = -p_kij * std::log(p_kij);
                population_entropy_k += e_kij;
            }
        }
        
        population_entropy[k] = population_entropy_k;
    }
}

arma::colvec Population::accumulated_proportional_fitness(
        const std::vector<Individual> & p
) {
    int N = p.size();
    arma::colvec accumulated_fitness(N);
    
    accumulated_fitness[0] = p[0].fitness_scaled;
    for (int n = 1; n < N; n++)
    {
        accumulated_fitness[n] = accumulated_fitness[n - 1] + p[n].fitness_scaled;
    }
    
    return accumulated_fitness / accumulated_fitness[N - 1];
}

//// Individual
// Constructors
Individual::Individual() 
{
    n_mutations = 0;
    
    age = 0.0;
    fitness = -HUGE_VAL;
    fitness_scaled = -HUGE_VAL;
}

Individual::Individual(
    const int & _n_genome, 
    const std::vector<int> & _n_chromosomes,
    Fitness & _fitness, 
    RV & _random_variate
) {
    n_mutations = 0;
    age = 0.0;
    
    chromosomes = std::vector<arma::colvec>(_n_genome);
    parameters = arma::colvec(_n_genome);
    for (int n = 0; n < _n_genome; n++)
    {
        const int & N = _n_chromosomes[n];
        
        arma::colvec chromosome = generate_random_chromosome(
            N, 
            _random_variate
        );
        
        chromosomes[n] = chromosome;
        
        const double & l = _fitness.lower[n];
        const double & u = _fitness.upper[n];
        
        parameters[n] = decode_chromosome(
            chromosome, 
            l, 
            u, 
            N
        );
    }
    
    fitness = _fitness.calculate_fitness(
        parameters
    );
    
    fitness_scaled = _fitness.calculate_fitness_scaled(
        fitness, 
        age, 
        n_mutations
    );
}

arma::colvec Individual::generate_random_chromosome(
        const int & N, 
        RV & random_variate
) {
    arma::colvec c(N);
    for (int n = 0; n < N; n++)
    {
        double u = random_variate.generate_uniform_real();
        if (u < 0.5) {
            c[n] = 0.0;
        }
        else 
        {
            c[n] = 1.0;
        }
    }
    
    return c;
}

double Individual::decode_chromosome(
        const arma::colvec & x, 
        const double & l, 
        const double & u, 
        const int & N
) {
    double sum = 0.0;
    for (int n = 0; n < N; n++) 
    {
        if (x[n] > 0.0) 
        {
            double bit = std::pow(2.0, n);
            sum += bit;
        }
    }
    
    const double & s = std::pow(2.0, N) - 1.0;
    const double & d = l + ((u - l) / s) * sum;
    
    return d;
}

void Individual::decode_individual(
        Fitness & fitness,
        const std::vector<int> & n_chromosomes,
        const int & n_genome
) {
    arma::colvec pars(n_genome);
    for (int n = 0; n < n_genome; n++) 
    {
        pars[n] = decode_chromosome(
            chromosomes[n], 
            fitness.lower[n], 
            fitness.upper[n],
            n_chromosomes[n]
        ); 
    }
    
    parameters = pars;
}
