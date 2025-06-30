#pragma once
#include <cstddef>

struct RunStats
{
    // --- runtimes ---
    double sat_time_total   = 0.0;   // sat operator runtime
    double other_time_total = 0.0;   // other operator runtime
    double overhead_total   = 0.0;   // framework overhead
    double validation_time_total = 0.0;
    double sat_ratio_ops    = 0.0;
    double sat_ratio        = 0.0;

    // --- count of calls/iterations ---
    std::size_t outer_iterations = 0; // outer (LNS) iteration
    std::size_t other_iters      = 0;
    std::size_t sat_iters        = 0; // how many times was SAT chosen as the operator
    std::size_t sat_calls        = 0; // how many times the SAT actually ran
    std::size_t sat_ok           = 0; // successful solution from SAT
    std::size_t sat_fail         = 0;

    // --- solution quality ---
    int initial_soc      = -1;
    int final_soc        = -1;
    int failed_iterations = 0;
    int sat_repairs         = 0;
    int sat_no_improve     = 0;

    // --- wall-time ---
    double wall_runtime = 0.0;
};