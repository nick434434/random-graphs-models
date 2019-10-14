//
// Created by nick434434 on 14-10-19.
//

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <vector>
#include "prob_stuff.h"

using std::vector;

static boost::mt19937 randGen(15); //held constant for repeatability
static double tau = uniform_distribution(generator) + 2 + 1e-10;
static boost::math::pareto_distribution<> dist(tau, tau);
static boost::random::uniform_real_distribution<> uniformReal(1.0,10.0); //this range can be adjusted to effect values
static boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > generator_uniform(randGen, uniformReal);


vector<long> pareto_vec(long n, long max) {
    vector<long> complements(n);
    for (long i = 0; i < n; i++)
        complements[i] = (long)((boost::math::cdf(boost::math::complement(dist, generator_uniform())) - 0.1) * 10 * max / 9);
    return complements;
}


vector<double> pareto_vec(long n) {
    vector<double> complements(n);
    for (long i = 0; i < n; i++)
        complements[i] = boost::math::cdf(boost::math::complement(dist, generator_uniform()));
    return complements;
}


std::function<long()> get_pareto_generator(long min, long max) {
    return [min, max]() {
        return (long)((boost::math::cdf(boost::math::complement(dist, generator_uniform())) - 0.1) * 10 * (max - min) / 9) + min;
    };
}
