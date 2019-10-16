//
// Created by nick434434 on 14-10-19.
//

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <vector>
#include <iostream>
#include "prob_stuff.h"

using std::vector;

static boost::mt19937 randGen(15); //held constant for repeatability
static double tau = uniform_distribution(generator) + 2 + 1e-10;
static boost::math::pareto_distribution<> dist(tau, tau);
static boost::random::uniform_real_distribution<> uniformReal(1.0,10.0); //this range can be adjusted to effect values
static boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > generator_uniform(randGen, uniformReal);


vector<long> pareto_vec(long n, long min, long max) {
    printf("Tau = %f\n", tau);
    vector<double> complements(n);
    for (long i = 0; i < n; i++)
        complements[i] = boost::math::cdf(boost::math::complement(dist, generator_uniform()));

    double actual_min = *std::min_element(complements.begin(), complements.end());
    double actual_max = *std::max_element(complements.begin(), complements.end());


    vector<long> counts(50);
    for (long i = 0; i < n; ++i)
        try {
            counts[(long)(std::trunc((complements[i] - actual_min) / (actual_max - actual_min) * 100)) / 2]++;
        } catch(std::out_of_range& e) {
            printf("%s\n", e.what());
        }

    for (long i = 0; i < counts.size(); ++i)
        std::cout << i << (i < 10 ? "  | " : " | ") << counts[i] << std::endl;
    std::cout << std::endl;

    vector<long> results;
    std::cout << "Maximum: " << actual_max << std::endl;
    // double step = (actual_max - actual_min) / (max - min + 1);
    long power = (long)std::trunc(log10(max - min + 1) + 1);
    power = std::pow(10, power);
    double cf = power / max;
    std::transform(complements.begin(), complements.end(), std::back_inserter(results), [actual_min, actual_max, min, max](double a) {
        long res = (long)(std::trunc((a - actual_min) / (actual_max - actual_min) * (max - min)) + min);
        //if (res >= max)
        //    std::cout << "FUCKERY: " << a << " became " << res << std::endl;
        return res;
    });
    return results;
}


vector<double> pareto_vec(long n) {
    printf("Tau = %f\n", tau);
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

