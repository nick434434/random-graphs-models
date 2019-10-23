//
// Created by nick434434 on 14-10-19.
//
#pragma once

#ifndef NUMERICALASSIGNMENTPS2_PROB_STUFF_H
#define NUMERICALASSIGNMENTPS2_PROB_STUFF_H

#include <random>
#include <cassert>
#include <vector>
#include <functional>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/math/common_factor.hpp>

using std::vector;


static std::random_device generator;
static std::uniform_real_distribution<double> uniform_distribution(0, 1);
auto edge_decider = [](double prob) {
    assert(prob > 0 && prob < 1);
    return uniform_distribution(generator) < prob;
};

static boost::mt19937 randGen(15); //held constant for repeatability
static double tau = uniform_distribution(generator) + 2 + 1e-10;
static boost::math::pareto_distribution<> dist(tau, tau);
static boost::random::uniform_real_distribution<> uniformReal(1.0,10.0); //this range can be adjusted to effect values
static boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > generator_uniform(randGen, uniformReal);


vector<long> pareto_vec(long n, long min, long max);

vector<double> pareto_vec(long n);

std::function<long()> get_pareto_generator(long min, long max);

vector<long> generate_n_pareto(long min, long max, long n);

vector<long> alternative_n_pareto(long min, long max, long n);

#endif //NUMERICALASSIGNMENTPS2_PROB_STUFF_H
