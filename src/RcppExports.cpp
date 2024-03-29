// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cbga_proportional
Rcpp::List cbga_proportional(Rcpp::Function& f, const arma::colvec& lower, const arma::colvec& upper, const int& n_lifetimes, const std::vector<int>& n_population, const std::vector<int> n_chromosomes, const double& pi_mutation, const double& pi_recombination, const int& max_lifespan, const double& s_age, const double& s_mutation, const int& breeding_rotation_rate, const arma::colvec& breeding_rotation, const int& trace);
RcppExport SEXP _CBGA_cbga_proportional(SEXP fSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP n_lifetimesSEXP, SEXP n_populationSEXP, SEXP n_chromosomesSEXP, SEXP pi_mutationSEXP, SEXP pi_recombinationSEXP, SEXP max_lifespanSEXP, SEXP s_ageSEXP, SEXP s_mutationSEXP, SEXP breeding_rotation_rateSEXP, SEXP breeding_rotationSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_lifetimes(n_lifetimesSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type n_population(n_populationSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type n_chromosomes(n_chromosomesSEXP);
    Rcpp::traits::input_parameter< const double& >::type pi_mutation(pi_mutationSEXP);
    Rcpp::traits::input_parameter< const double& >::type pi_recombination(pi_recombinationSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_lifespan(max_lifespanSEXP);
    Rcpp::traits::input_parameter< const double& >::type s_age(s_ageSEXP);
    Rcpp::traits::input_parameter< const double& >::type s_mutation(s_mutationSEXP);
    Rcpp::traits::input_parameter< const int& >::type breeding_rotation_rate(breeding_rotation_rateSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type breeding_rotation(breeding_rotationSEXP);
    Rcpp::traits::input_parameter< const int& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(cbga_proportional(f, lower, upper, n_lifetimes, n_population, n_chromosomes, pi_mutation, pi_recombination, max_lifespan, s_age, s_mutation, breeding_rotation_rate, breeding_rotation, trace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CBGA_cbga_proportional", (DL_FUNC) &_CBGA_cbga_proportional, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_CBGA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
