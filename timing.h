//
// Created by nick434434 on 23/10/2019.
//
#pragma once

#ifndef RANDOM_GRAPHS_MODELS_TIMING_H
#define RANDOM_GRAPHS_MODELS_TIMING_H

#include <chrono>
#include <iostream>

namespace timing {
    struct Time {
        long microsec;
        Time(): microsec(0) {}
        explicit Time(long val): microsec(val) {}
        double as_seconds() { return (double)microsec / 1e6; }
        double as_milliseconds() { return (double)microsec / 1e3; }
        double as_minutes() { return (double)microsec / 6e7; }
    };

    static auto t0 = std::chrono::high_resolution_clock::now();

    static auto t1 = std::chrono::high_resolution_clock::now();

    void start_clock();

    Time reset_clock(std::string& message, std::ostream& out);

    Time check_clock();

    void start_local_clock();

    Time reset_local_clock(std::string& message, std::ostream& out);

    Time check_local_clock();
}

std::ostream& operator<<(std::ostream& out, timing::Time t);

#endif //RANDOM_GRAPHS_MODELS_TIMING_H
