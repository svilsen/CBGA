#' @title Controlled Breeding Genetic Algorithm Control
#' 
#' @description \code{cbga.control} creates a list of all parameters needed to run the genetic algorithm. 
#' 
#' @param n_lifetime The number of iterations of the genetic algorithm.
#' @param n_population A vector setting number of individuals in the active populations (the size of the inactive population is the sum this vector).
#' @param n_chromosomes A vector setting the length of the bitstrings used for each individual (the length of the bitstrings relates to the precision of the optimisation).
#' @param pi_mutation The probability of mutation.
#' @param pi_recombination The probability of recombination.
#' @param s_age A scaling factor for the effect of age on the fitness (> 0).
#' @param s_mutation A scaling factor for the effect of the number of mutations on the fitness (> 0).
#' @param max_lifespan The maximum lifespan of an individual (0 < \code{max_lifespan} < \code{n_lifetime}).
#' @param selection_type A character string setting the type of selection used in litter, and active population selection (default is \code{"proportional"}).
#' @param trace An integer between 0 and \code{n_lifetime}. If zero then trace is not shown. If larger than 0, then a trace is shown every \code{trace} iterations.
#' 
#' @return A list of parameters used for the \link{cbga}-function.
#' @examples
#' cbga.control()
#' 
#' @export
cbga.control <- function(
    n_lifetime = 1000, 
    n_population = c(10, 10, 10, 10, 10, 10), 
    n_chromosomes = NULL,
    pi_mutation = NULL, 
    pi_recombination = NULL,
    s_age = 4, 
    s_mutation = 1, 
    max_lifespan = 10,
    selection_type = "proportional",
    trace = 0
) {
    if (is.null(n_lifetime) || !is.numeric(n_lifetime)) {
        n_lifetime <- 100
    }
    
    if (is.null(n_population) || !is.numeric(n_population)) {
        n_population <- c(10, 10, 10, 10, 10, 10)
    }
    
    if (is.null(n_chromosomes)|| !is.numeric(n_chromosomes)) {
        n_chromosomes <- 20
    }
    
    if (is.null(s_age) || !is.numeric(s_age) || (s_age < 0)) {
        s_age <- 4
    }
    
    if (is.null(s_mutation) || !is.numeric(s_mutation) || (s_mutation < 0)) {
        s_mutation <- 1
    }
    
    if (is.null(max_lifespan)) {
        max_lifespan <- min(5, round(n_lifetimes / 10))
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
    
    if (is.null(selection_type) || (selection_type %in% c("p", "prop", "proportional"))) {
        selection_type <- "proportional"
    }
    
    control <- list(n_lifetime = n_lifetime, n_population = n_population, n_chromosomes = n_chromosomes, 
                    pi_mutation = pi_mutation, pi_recombination = pi_recombination, 
                    max_lifespan = max_lifespan, s_age = s_age, s_mutation = s_mutation,  
                    selection_type = selection_type, trace = trace)
    return(control)
}

#' @title Controlled Breeding Genetic Algorithm
#' 
#' @description \code{cbga} is a genetic algorithm used to fit box constrained univariate functions. The genetic algorithm utilises controlled breeding for a more efficient estimation.
#' 
#' @param f The function to be optimised by the algorithm.
#' @param lower A vector of lower bounds.
#' @param upper A vector of upper bounds.
#' @param bp A list of arguments passed to the \link{breeding_protocol}-function (see \code{?breeding_protocol} for more information).
#' @param ... Additional control parameters used in the genetic algorithm. See \link{cbga.control} for more information.
#' 
#' @return The fittest individual found by the genetic algorithm, arranged in the list containing the following:
#' 
#' @examples 
#' \dontrun{
#' f <- function(x) exp(-sum(x)^2)
#' lower <- rep(-1, 3)
#' upper <- rep(1, 3)
#' 
#' cbga(f, lower, upper)
#' }
#' 
#' @export
cbga <- function(
    f, 
    lower, 
    upper, 
    bp = list(), 
    ...
) {
    ## Control parameters
    control <- list(...)
    control <- do.call(cbga.control, control)
    if (length(lower) != length(upper)) {
        stop("The 'lower' and 'upper' bounds have to be the same length.")
    }
    
    if (length(lower) != length(control$n_chromosomes)) {
        if (length(control$n_chromosomes) == 1) {
            control$n_chromosomes <- rep(control$n_chromosomes, length(lower))
        }
        else {
            stop("'n_chromosomes' needs to have the same length as the bounds.")
        }
    }
    
    if (is.null(control$pi_mutation)) {
        control$pi_mutation <- 1.0 / (sum(control$n_chromosomes))
    }
    
    if (is.null(control$pi_recombination)) {
        control$pi_recombination <- 1.0 / (sum(control$n_chromosomes))
    }
    
    ## Breeding parameters
    if (is.null(bp)) {
        bp <- list()
    }
    
    bp$n_population <- control$n_population
    bp$n_lifetime <- control$n_lifetime
    bp <- do.call(breeding_protocol, bp)
    
    ## CBGA
    if (control$selection_type == "proportional") {
        res <- cbga_proportional(
            f = f, 
            lower = matrix(lower, ncol = 1), 
            upper = matrix(upper, ncol = 1), 
            n_lifetimes = control$n_lifetime, 
            n_population = control$n_population, 
            n_chromosomes = matrix(control$n_chromosomes, ncol = 1), 
            pi_mutation = control$pi_mutation, 
            pi_recombination = control$pi_recombination, 
            max_lifespan = control$max_lifespan,
            s_age = control$s_age, 
            s_mutation = control$s_mutation, 
            breeding_rotation_rate = bp$breeding_rotation_rate,
            breeding_rotation = bp$breeding_rotation,
            trace = control$trace
        )
    }
    else {
        stop("Only implemented selection type is proportional.")
    }
    
    ## 
    class(res) <- "cbga"
    return(res)
}

#' @export
coef.cbga <- function(
    object, 
    ...
) {
    object$parameters
}


