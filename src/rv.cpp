#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]

#include <ctime>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include "rv.hpp"

rv::rv() : 
    rng(std::time(0)), uniform_real(), generate_uniform_real(rng, uniform_real) { }
