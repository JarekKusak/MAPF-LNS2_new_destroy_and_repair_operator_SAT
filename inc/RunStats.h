#pragma once
#include <cstddef>

struct RunStats
{
    // --- kumulované časy ---
    double sat_time_total   = 0.0;   // celkový čas strávený v SAT solveru
    double other_time_total = 0.0;   // re-plan (PP/CBS/…) + init-solution
    double overhead_total   = 0.0;   // rámcový “overhead”
    double sat_ratio_ops    = 0.0;
    double sat_ratio        = 0.0;

    // --- počty volání / iterací ---
    std::size_t outer_iterations = 0; // vnější (LNS) iterace
    std::size_t sat_iters        = 0; // kolikrát byl jako operátor zvolen SAT
    std::size_t sat_calls        = 0; // kolikrát se SAT skutečně spustil
    std::size_t sat_ok           = 0; // úspěšné řešení od SAT
    std::size_t sat_fail         = 0;

    // --- kvalita výsledku ---
    int initial_soc      = -1;
    int final_soc        = -1;
    int failed_iterations = 0;

    // --- celkový wall-time ---
    double wall_runtime = 0.0;
};