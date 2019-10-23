//
// Created by nick434434 on 23/10/2019.
//

#include "timing.h"

void start_clock() {
    t0 = std::chrono::high_resolution_clock::now();
}

long reset_clock() {
    auto diff = std::chrono::high_resolution_clock::now() - t0;
    t0 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(diff).count();
}

long check_clock() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0).count();
}
