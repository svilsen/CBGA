#include <RcppArmadillo.h>

#include "rv.hpp"
#include "individual.hpp"

//// Fitness
// Constructor
Fitness::Fitness(Rcpp::Function & _f, const arma::colvec & _lower, const arma::colvec & _upper, 
                 const double & _s_age, const double & _pi_mutation, const double & _e_mutation, 
                 const double & _s_mutation) : 
    f(_f), lower(_lower), upper(_upper), s_age(_s_age), 
    pi_mutation(_pi_mutation), e_mutation(_e_mutation), s_mutation(_s_mutation) { }

// Functions
double Fitness::age_scaling(const double & age)
{
    double age_s = s_age * age * (1.0 - age);
    if (age_s > 1.0) 
    {
        age_s = 1.0;
    }
    
    return age_s;
}

double Fitness::mutation_scaling(const int & n_mutations)
{
    double v = std::pow(pi_mutation * (1.0 - pi_mutation), 0.5);
    double r_standard = (n_mutations - e_mutation) / v;
    
    return std::exp(-s_mutation * std::abs(r_standard));
}

double Fitness::decode_chromosome(const arma::colvec & x, const double & l, const double & u, const int & N) 
{
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

void Fitness::decode_individual(Individual & individual, const std::vector<int> & n_chromosomes, const int & n_gnome) 
{
    arma::colvec pars(n_gnome);
    for (int n = 0; n < n_gnome; n++) 
    {
        pars[n] = decode_chromosome(individual.chromosomes[n], lower[n], upper[n], n_chromosomes[n]); 
    }
    
    individual.parameters = pars;
}

double Fitness::calculate_fitness(const arma::colvec & paramters, const double & age, const int & n_mutations) 
{
    double age_s = age_scaling(age);
    double mutation_s = mutation_scaling(n_mutations, e_mutation, pi_mutation);
    double fitness = f(parameters);
    
    return age_s * mutation_s * fitness;
}


//// Population 
// Constructor
population::population(const int & _n_population, const int & _n_genome, const std::vector<int> _n_chromosome, 
                       const Fitness & _fitness, const RV & _random_variate) :
    n_population(_n_population), n_genome(_n_genome), n_chromosome(_n_chromosome)
{
    std::vector<Individual> p_active(n_genome);
    for (int n = 0; n < n_genome; n++) 
    {
        Individual individual(n_genome, n_chromosomes, _fitness, _random_variate);
        p_active[n] = individual;
    }
    
    update_inbreeding_coefficient();
}

void population::update_inbreeding_coefficient()
{
    
}

//// Individual
// Constrctors
Individual::Individual() 
{
    n_mutations = 0;
    
    age = 0.0;
    fitness = -HUGE_VAL;
}

Individual::Individual(const int & _n_genome, const std::vector<int> & _n_chromosomes, 
                       const Fitness & _fitness, const RV & _random_variate) 
{
    n_mutations = 0;
    age = 0.0;
    
    chromosomes = std::vector<arma::colvec>(_n_genome);
    parameters = arma::colvec(_n_gnome);
    for (int n = 0; n < _n_gnome; n++)
    {
        int & N = _n_chromosomes[n];
        arma::colvec chromosome = generate_random_chromosome(N, _random_variate);
        chromosomes[n] = chromosome;
        
        double & l = _fitness.lower[n];
        double & u = _fitness.upper[n];
        parameters[n] = _fitness.decode_chromosome(chromosome, l, u, N);
    }
    
    fitness = _fitness.calculate_fitness(parameters, age, n_mutations);
}

arma::colvec Individual::generate_random_chromosome(const int & N, const RV & random_variate) 
{
    arma::colvec c = arma::colvec::zeros(N);
    for (int n = 0; n < N; n++)
    {
        double u = random_variate.generate_uniform_real();
        if (u > 0.5) {
            c[n] = 1.0;
        }
    }
    
    return c;
}
