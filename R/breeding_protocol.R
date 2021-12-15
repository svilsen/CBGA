#' @title Breeding protocols
#' 
#' @description \code{breeding_protocol} creates a list of all parameters needed to control the breeding schedule of the \link{cbga} algorithm.
#' 
#' @param n_lifetime The number of iterations of the genetic algorithm.
#' @param n_population A vector setting number of individuals in the active populations (the size of the inactive population is the sum this vector).
#' @param breeding_rotation_rate The rate at which the breeding rotation schedule is applied in the protocol. 
#' @param breeding_rotation A vector indicating the rotation of the bulls in the population. 
#' 
#' @return A list containing the 'breeding_rotation_rate' and 'breeding_rotation'.
#' @export 
breeding_protocol <- function(n_population, n_lifetime, breeding_rotation_rate = 2, breeding_rotation = NULL) {
    ## Rotation rate
    if (is.null(breeding_rotation_rate)) {
        breeding_rotation_rate <- 2
    }
    
    if (breeding_rotation_rate < 1) {
        warning("'breeding_rotation_rate' was set smaller than 1, setting 'breeding_rotation_rate' to 2.")
        breeding_rotation_rate <- 2
    }
    
    if (breeding_rotation_rate > n_lifetime) {
        warning("'breeding_rotation_rate' was set larger than 'n_lifetime', setting 'breeding_rotation_rate' to 2.")
        breeding_rotation_rate <- 2
    }
    
    ## Rotation schedule
    K <- length(n_population)
    
    if (is.null(breeding_rotation) || breeding_rotation == "random") {
        breeding_rotation <- sample(seq_len(K), K, replace = FALSE)
    } 
    
    if (!is.numeric(breeding_rotation)) {
        warning("'breeding_rotation' is not a numeric vector, randomising 'breeding_rotation'.")
        breeding_rotation <- sample(seq_len(K), K, replace = FALSE)
    } 
    
    if (all(!(breeding_rotation %in% seq_len(K)))) {
        stop("'breeding_rotation' points to invalid elements, it has to values in the set {1, 2, ..., K}.")
    }
    
    breeding_rotation <- matrix(breeding_rotation - 1, ncol = 1)
    
    ## 
    return(list(breeding_rotation_rate = breeding_rotation_rate, breeding_rotation = breeding_rotation))
}