//
// Created by Jaroslav Kusák on 09.06.2025.
//

#pragma once
#include <iostream>
#include <atomic>

namespace satlog {
    // runtime-flag – set from CLI
    inline std::atomic_bool debug_enabled {true};   // default: ON (for script)

    /* macros */
    #define SAT_STAT(msg)  do { std::cout << "[STAT] " << msg << std::endl; } while(0)
    #define SAT_DBG(msg)   do { if (satlog::debug_enabled) \
                                       std::cout << "[DBG]  " << msg << std::endl; } while(0)
}