//
// Created by nick434434 on 14-10-19.
//
#pragma once

#ifndef NUMERICALASSIGNMENTPS2_ERDOS_RENYI_H
#define NUMERICALASSIGNMENTPS2_ERDOS_RENYI_H

#include <vector>

using std::pair;
using std::vector;

void print_graph(vector<vector<int>> G);

vector<vector<int>> generate_ER(int n, double lambd);

long count_edges(vector<vector<int>> graph);

long count_triangles(vector<vector<int>> graph);

pair<vector<double>, vector<long>> generate_for_single_n(int N, double lambda_start, double lambda_step, int num_iters);

pair<vector<int>, vector<long>> generate_for_single_lambda(double lmbd, int N_start, int N_step, int num_iters);

vector<double> generate_poisson_function(double lmbd, long start, long end, long scale);

#endif //NUMERICALASSIGNMENTPS2_ERDOS_RENYI_H
