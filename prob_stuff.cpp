//
// Created by nick434434 on 14-10-19.
//

#include <iostream>
#include "prob_stuff.h"

using std::vector;


vector<long> pareto_vec(long n, long min, long max) {
    printf("Tau = %f\n", tau);
    vector<double> complements(n);
    double val;
    double eps = 1e-10;
    for (long i = 0; i < n; i++) {
        while (abs((val = boost::math::cdf(boost::math::complement(dist, generator_uniform()))) - 1) < eps);
        complements[i] = val;
    }

    double actual_min = *std::min_element(complements.begin(), complements.end());
    double actual_max = *std::max_element(complements.begin(), complements.end());
    std::cout << "Min max: " << actual_min << "; " << actual_max << std::endl;

    vector<long> results;
    std::transform(complements.begin(), complements.end(), std::back_inserter(results), [actual_min, actual_max, min, max](double a) {
        return (long)std::trunc((a - actual_min) / (actual_max - actual_min) * (max - min)) + min;
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
        double val;
        while (abs((val = boost::math::cdf(boost::math::complement(dist, generator_uniform()))) - 1) < 1e-10);
        return (long)((val - 0.1) * 10 * (max - min) / 9) + min;
    };
}


std::function<double()> get_pareto_generator_std() {
    return []() {
        double val;
        while (abs((val = boost::math::cdf(boost::math::complement(dist, generator_uniform()))) - 1) < 1e-10);
        return val;
    };
}


vector<long> generate_n_pareto(long min, long max, long n) {
    vector<long> res;
    vector<double> raw(n);
    auto gen = get_pareto_generator_std();
    for (long i = 0; i < n; ++i)
        raw[i] = gen();

    double actual_min = *std::min_element(raw.begin(), raw.end());
    double actual_max = *std::max_element(raw.begin(), raw.end());

    std::transform(raw.begin(), raw.end(), std::back_inserter(res), [actual_min, actual_max, min, max](double a) {
        return (long)std::trunc((a - actual_min) / (actual_max - actual_min) * (max - min)) + min;
    });
    return res;
}


vector<long> alternative_n_pareto(long n) {
    vector<long> res(n);

    for (long i = 0; i < n; ++i)
        res[i] = (long)std::pow(1. / uniform_distribution(generator), 1. / (tau - 1));

    return res;
}


vector<long> alternative_n_pareto_truncated(long n, double gamma, bool beta_loglog) {
    vector<double> res(n);
    double beta = 1. / std::pow(beta_loglog ? log(log(n)) : log(n), gamma);
    long n_beta = (long)std::pow(n, beta);

    printf("beta = %1.4f, n = %ld, n_beta = %ld\n", beta, n, n_beta);

    for (long i = 0; i < n; ++i)
        res[i] = uniform_distribution(generator);

    vector<long> transformed;
    std::transform(res.begin(), res.end(), std::back_inserter(transformed), [n_beta](double a) {
        long res = (long)std::pow(1 / a, 1 / (tau - 1));
        while (res > n_beta)
            res = (long)std::pow(1 / uniform_distribution(generator), 1 / (tau - 1));
        return res;
    });

    return transformed;
}

