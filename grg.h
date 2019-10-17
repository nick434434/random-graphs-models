//
// Created by nick434434 on 14-10-19.
//
#pragma once

#ifndef NUMERICALASSIGNMENTPS2_GRG_H
#define NUMERICALASSIGNMENTPS2_GRG_H

#include <vector>
#include <functional>


vector<vector<int>> generate_GRG(long n, long max_weight, double tau_subs = 0);

vector<pair<long, long>> generate_GRG_edge_pairs(long n, const vector<long>& w);

vector<pair<long, long>> generate_GRG_edge_pairs(long n, long max_weight, double tau_subs = 0);

vector<vector<pair<long, long>>> generate_m_GRG_edge_pairs_same_weights(long n, long m, long max_weight, vector<long>& weights, double tau_subs = 0);

vector<vector<pair<long, long>>> generate_m_GRG_edge_pairs(long n, long m, long max_weight, double tau_subs = 0);

vector<long> count_degrees(const vector<pair<long, long>>& G);

#endif //NUMERICALASSIGNMENTPS2_GRG_H
