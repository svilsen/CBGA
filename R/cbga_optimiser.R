#' @title Controlled Breeding Genetic Algorithm Control
#' 
#' @description \code{cbga.control} creates a list of all parameters needed to run the genetic algorithm. 
#' 
#' @param n_lifetimes The number of iterations of the genetic algorithm.
#' @param n_population The number of individuals in both the active and inactive populations.
#' @param n_chromosomes A vector assigning the length of the bitstrings (related to the precision of the optimisation).
#' @param pi_mutation The probability of mutation.
#' @param pi_recombination The probability of recombination.
#' @param s_age A scaling factor for the effect of age on the fitness (> 0).
#' @param s_mutation A scaling factor for the effect of the number of mutations on the fitness (> 0).
#' @param trace TRUE/FALSE: Should an internal trace of the genetic algorithm be shown? 
#' 
#' @return A list of parameters used for the \link{cbga}-function.
#' @export
cbga.control <- function(n_lifetimes = 100, n_population = 100, n_chromosomes = NULL, 
                         pi_mutation = NULL, pi_recombination = NULL, 
                         s_age = 4, s_mutation = 1, 
                         trace = TRUE) {
    if (is.null(n_lifetimes) || !is.numeric(n_lifetimes)) {
        s_age <- 100
    }
    if (is.null(n_population) || !is.numeric(n_population)) {
        n_population <- 100
    }
    if (is.null(s_age) || !is.numeric(s_age) || (s_age < 0)) {
        s_age <- 4
    }
    if (is.null(s_mutation) || !is.numeric(s_mutation) || (s_mutation < 0)) {
        s_mutation <- 1
    }
    if (is.null(trace) || !is.logical(trace)) {
        trace <- FALSE
    }
    
    list(n_lifetimes = n_lifetimes, n_population = n_population, n_chromosomes = n_chromosomes, 
         pi_mutation = pi_mutation, pi_recombination = pi_recombination, 
         s_age = s_age, s_mutation = s_mutation, 
         trace = trace)
}

#' @title Controlled Breeding Genetic Algorithm
#' 
#' @description \code{cbga} is a genetic algorithm used to fit box contrained univarite functions. The genetic algorithm utilises controlled breeding for a more efficient estimation.
#' 
#' @param f The function to be optimised.
#' @param lower A vector of lower bounds.
#' @param upper A vector of upper bounds.
#' @param bp_type The breeding protocol used in the genetic algorithm. See details for more information.
#' @param ga_pars A list of additional parameters used in the genetic algorithm. See \link{cbga.control} for more information.
#' 
#' @details The following breeding protocols...
#' 
#' @return The fittest individual found by the genetic algorithm.
#' @export
cbga <- function(f, lower, upper, bp_type = "proportional", ga_pars = list()) {
    ga_pars <- do.call(cbga.control, ga_pars)
    if (length(lower) != length(upper)) {
        stop("The 'lower' and 'upper' bounds have to be the same length.")
    }
    
    if (is.null(ga_pars$n_chromosomes)) {
        ga_pars$n_chromosomes <- rep(10, length(lower))
    }
    if (length(lower) != length(ga_pars$n_chromosomes)) {
        stop("'n_chromosomes' needs to have the same length as the bounds. Set using the 'ga_pars' argument.")
    }
    
    if (is.null(ga_pars$pi_mutation)) {
        ga_pars$pi_mutation <- 1.0 / (sum(ga_pars$n_chromosomes))
    }
    if (is.null(ga_pars$pi_recombination)) {
        ga_pars$pi_recombination <- 1.0 / (sum(ga_pars$n_chromosomes))
    }
    
    if (tolower(bp_type) %in% c("p", "prop", "proportional")) {
        res <- cbga_proportional(f = f, 
                                 lower = lower, 
                                 upper = upper, 
                                 n_lifetimes = ga_pars$n_lifetimes, 
                                 n_population = ga_pars$n_population, 
                                 n_chromosomes = ga_pars$n_chromosomes, 
                                 pi_mutation = ga_pars$pi_mutation, 
                                 pi_recombination = ga_pars$pi_recombination, 
                                 s_age = ga_pars$s_age, 
                                 s_mutation = ga_pars$s_mutation, 
                                 trace = ga_pars$trace)
    }
    else {
        stop(paste0("The supplied breeding protocol '", bp_type ,"' is not a valid option. See details of '?cbga' for more information."))
    }

    class(res) <- "cbga"
    return(res)
}

#' @export
coef.cbga <- function(object, ...) {
    object$parameters
}
    
    
    