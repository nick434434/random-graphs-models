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

using std::vector;

static std::random_device generator;
static std::uniform_real_distribution<double> uniform_distribution(0, 1);
auto edge_decider = [](double prob) {
    assert(prob > 0 && prob < 1);
    return uniform_distribution(generator) < prob;
};

vector<long> pareto_vec(long n, long max);

vector<double> pareto_vec(long n);

std::function<long()> get_pareto_generator(long min, long max);

#endif //NUMERICALASSIGNMENTPS2_PROB_STUFF_H
