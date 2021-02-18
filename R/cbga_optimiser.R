#' @title Controlled Breeding Genetic Algorithm Control
#' 
#' @description \code{cbga.control} creates a list of all parameters needed to run the genetic algorithm. 
#' 
#' @param n_lifetime The number of iterations of the genetic algorithm.
#' @param n_population The number of individuals in both the active and inactive populations.
#' @param n_chromosomes A vector setting the length of the bitstrings used for each individual (the length of the bitstrings relates to the precision of the optimisation).
#' @param pi_mutation The probability of mutation.
#' @param pi_recombination The probability of recombination.
#' @param s_age A scaling factor for the effect of age on the fitness (> 0).
#' @param s_mutation A scaling factor for the effect of the number of mutations on the fitness (> 0).
#' @param max_lifespan The maximum lifespan of an individual (0 < \code{max_lifespan} < \code{n_lifetime}).
#' @param trace An integer between 0 and \code{n_lifetime}. If zero then trace is not shown. If larger than 0, then a trace is shown every \code{trace} iterations.
#' 
#' @return A list of parameters used for the \link{cbga}-function.
#' @examples
#' cbga.control()
#' 
#' @export
cbga.control <- function(n_lifetime = 100, n_population = 100, n_chromosomes = NULL, 
                         pi_mutation = NULL, pi_recombination = NULL, 
                         s_age = 4, s_mutation = 1, max_lifespan = 10,
                         trace = 0) {
    if (is.null(n_lifetime) || !is.numeric(n_lifetime)) {
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
    
    if (is.null(max_lifespan)) {
        max_lifespan <- round(n_lifetimes / 10)
    }
    if (max_lifespan < 1) {
        warning("'max_lifespan' has to be larger than 0, and is, therefore, set to 1.")
        max_lifespan <- 1
    }
    if (max_lifespan >= n_lifetime) {
        warning("'max_lifespan' has to be smaller than 'n_lifetime', and is, therefore, set to 'n_lifetime' - 1.")
        max_lifespan <- n_lifetime - 1
    }
    
    if (is.null(trace)) {
        trace <- 0
    }
    
    
    control <- list(n_lifetime = n_lifetime, n_population = n_population, n_chromosomes = n_chromosomes, 
                    pi_mutation = pi_mutation, pi_recombination = pi_recombination, 
                    max_lifespan = max_lifespan, s_age = s_age, s_mutation = s_mutation, 
                    trace = trace)
    return(control)
}

#' @title Controlled Breeding Genetic Algorithm
#' 
#' @description \code{cbga} is a genetic algorithm used to fit box constrained univariate functions. The genetic algorithm utilises controlled breeding for a more efficient estimation.
#' 
#' @param f The function to be optimised by the algorithm.
#' @param lower A vector of lower bounds.
#' @param upper A vector of upper bounds.
#' @param bp A string, or \link{bp}-object, defining the breeding protocol used in the genetic algorithm. See details for more information on allowed pre-sets supplied as strings, and \link{set_bp} for more information on possible breeding protocols.
#' @param control A list of control parameters used in the genetic algorithm. See \link{cbga.control} for more information.
#' 
#' @details The following breeding protocols...
#' 
#' @return The fittest individual found by the genetic algorithm, arranged in the list containing the following:
#' 
#' 
#' @examples 
#' \dontrun{
#' f <- function(x) exp(-x^2)
#' lower <- -1
#' upper <- 1
#' 
#' cbga(f, lower, upper)
#' }
#' 
#' @export
cbga <- function(f, lower, upper, bp = "proportional", control = list()) {
    ga_pars <- do.call(cbga.control, control)
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
    
    if (tolower(bp) %in% c("p", "prop", "proportional")) {
        res <- cbga_proportional(f = f, 
                                 lower = matrix(lower, ncol = 1), 
                                 upper = matrix(upper, ncol = 1), 
                                 n_lifetime = ga_pars$n_lifetime, 
                                 n_population = ga_pars$n_population, 
                                 n_chromosomes = matrix(ga_pars$n_chromosomes, ncol = 1), 
                                 pi_mutation = ga_pars$pi_mutation, 
                                 pi_recombination = ga_pars$pi_recombination, 
                                 max_lifespan = ga_pars$max_lifespan,
                                 s_age = ga_pars$s_age, 
                                 s_mutation = ga_pars$s_mutation, 
                                 trace = ga_pars$trace)
    }
    else {
        stop(paste0("The supplied breeding protocol is not a valid. See details of '?cbga' and '?set_bp' for more information."))
    }
    
    class(res) <- "cbga"
    return(res)
}

#' @export
coef.cbga <- function(object, ...) {
    object$parameters
}


