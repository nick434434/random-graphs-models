//
// Created by nick434434 on 23/10/2019.
//
#pragma once

#ifndef RANDOM_GRAPHS_MODELS_TIMING_H
#define RANDOM_GRAPHS_MODELS_TIMING_H

#include <chrono>

static auto t0 = std::chrono::high_resolution_clock::now();

void start_clock();

long reset_clock();

long check_clock();

#endif //RANDOM_GRAPHS_MODELS_TIMING_H
