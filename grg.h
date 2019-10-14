//
// Created by nick434434 on 14-10-19.
//
#pragma once

#ifndef NUMERICALASSIGNMENTPS2_GRG_H
#define NUMERICALASSIGNMENTPS2_GRG_H

#include<vector>
#include <functional>

using std::vector;

vector<vector<int>> generate_GRG(long n, const std::function<long()>&);

vector<vector<int>> generate_GRG(long n, vector<long> weights);

#endif //NUMERICALASSIGNMENTPS2_GRG_H
