#ifndef rv_hpp
#define rv_hpp

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class rv 
{ 
public:
    boost::mt19937 rng;
    boost::random::uniform_01<> uniform_real;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > generate_uniform_real;
    
    rv();
};

#endif