//
// Created by nick434434 on 23/10/2019.
//

#include "timing.h"
#include <iostream>
#include <iomanip>

void timing::start_clock() {
    t0 = std::chrono::high_resolution_clock::now();
}

timing::Time timing::reset_clock(std::string& message, std::ostream& out) {
    if (message.length() > 40)
        message.resize(40);
    auto diff = std::chrono::high_resolution_clock::now() - t0;
    t0 = std::chrono::high_resolution_clock::now();
    Time tt = Time(std::chrono::duration_cast<std::chrono::microseconds>(diff).count());
    out << "Global clock reset due to \"" << std::setw(40) << std::left << message << "\" on " <<
        std::setw(10) << std::setprecision(2) << std::left << std::fixed << tt << " seconds" << std::endl;
    return tt;
}

timing::Time timing::check_clock() {
    return Time(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0).count());
}

void timing::start_local_clock() {
    t1 = std::chrono::high_resolution_clock::now();
}

timing::Time timing::reset_local_clock(std::string& message, std::ostream& out) {
    if (message.length() > 40)
        message.resize(40);
    auto diff = std::chrono::high_resolution_clock::now() - t1;
    t1 = std::chrono::high_resolution_clock::now();
    Time tt = Time(std::chrono::duration_cast<std::chrono::microseconds>(diff).count());
    out << "Local clock reset due to \"" << std::setw(40) << std::left << message << "\" on " <<
        std::setw(10) << std::setprecision(2) << std::left << std::fixed << tt << " seconds" << std::endl;
    return tt;
}

timing::Time timing::check_local_clock() {
    return Time(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t1).count());
}

std::ostream& operator<<(std::ostream& op_out, timing::Time t) {
    op_out << t.as_seconds(); return op_out;
}