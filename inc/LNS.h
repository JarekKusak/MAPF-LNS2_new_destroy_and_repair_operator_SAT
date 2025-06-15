#pragma once
#include "BasicLNS.h"
#include "InitLNS.h"

//pibt related
#include "simplegrid.h"
#include "pibt_agent.h"
#include "problem.h"
#include "mapf.h"
#include "pibt.h"
#include "pps.h"
#include "winpibt.h"

#include <random>

enum destroy_heuristic { RANDOMAGENTS, RANDOMWALK, INTERSECTION, DESTORY_COUNT, SAT, NONE };

// TODO: adaptively change the neighbor size, that is,
// increase it if no progress is made for a while
// decrease it if replanning fails to find any solutions for several times

#include <stdexcept>
#include <string>

class ValidationException : public std::runtime_error {
public:
    explicit ValidationException(const std::string& msg) : std::runtime_error(msg) {}
};


class LNS : public BasicLNS
{
public:
    vector<Agent> agents;
    double preprocessing_time = 0;
    double initial_solution_runtime = 0;
    int initial_sum_of_costs = -1;
    int sum_of_costs_lowerbound = -1;
    int sum_of_distances = -1;
    int restart_times = 0;

    LNS(const Instance& instance, double time_limit,
        const string & init_algo_name, const string & replan_algo_name, const string & destory_name,
        int neighbor_size, int num_of_iterations, bool init_lns, const string & init_destory_name, bool use_sipp,
        int screen, PIBTPPS_option pipp_option, const string& sat_heur_name, int sat_submap_side, int sat_prob_percent,
        const std::string& fallback_dest_name,
        const std::string& fallback_replan_name);
    ~LNS()
    {
        delete init_lns;
    }
    bool getInitialSolution();
    bool run();
    void validateSolution() const;
    void writeIterStatsToFile(const string & file_name) const;
    void writeResultToFile(const string & file_name) const;
    void writePathsToFile(const string & file_name) const;
    string getSolverName() const override { return "LNS(" + init_algo_name + ";" + replan_algo_name + ")"; }
private:
    InitLNS* init_lns = nullptr;
    string init_algo_name;
    string replan_algo_name;
    bool use_init_lns; // use LNS to find initial solutions
    destroy_heuristic destroy_strategy = RANDOMWALK;
    int num_of_iterations;
    string init_destory_name;
    PIBTPPS_option pipp_option;


    PathTable path_table; // 1. stores the paths of all agents in a time-space table;
    // 2. avoid making copies of this variable as much as possible.
    unordered_set<int> tabu_list; // used by randomwalk strategy
    list<int> intersections;

    std::set<std::pair<int,int>> ignored_agents_with_timestep;

    bool runEECBS();
    bool runCBS();
    bool runPP();
    bool runPIBT();
    bool runPPS();
    bool runWinPIBT();


    MAPF preparePIBTProblem(vector<int>& shuffled_agents);
    void updatePIBTResult(const PIBT_Agents& A, vector<int>& shuffled_agents);

    void chooseDestroyHeuristicbyALNS();

    bool generateNeighborByRandomWalk();
    bool generateNeighborByIntersection();

    int findMostDelayedAgent(); // the function is now returning the most problematic agent and his most problematic time step
    int findRandomAgent() const;
    void randomWalk(int agent_id, int start_location, int start_timestep,
                    set<int>& neighbor, int neighbor_size, int upperbound);

    bool generateNeighborBySAT(); // new destroy operator
    bool runSAT(); // new repair operator
    pair<vector<vector<int>>, vector<int>> getSubmapAndAgents(int agent_id, int submap_size, int agent_location, int timestep); // helper function for getting sub-map
    int last_selected_agent = -1;   // udržuje, na kom jsme skončili
    int countConflicts(const Agent& ag);
    int computeMaxDelay(const Agent& ag);
    double agentScore(const Agent& ag) const;
    pair<int,int> findBestAgentAndTime();
    const double W_DELAY_init   = 0.25;//4.0;
    const double W_CONFL_init   = 0.25;//2.0;
    const double W_STRETCH_init = 0.25;//1.0;
    const double W_REC_init     = 0.25;//0.5;
    void updateAllStats(int iter);
    int current_iter;
    void updateComponentWeights(int metric_index, double delta);
    std::vector<double> component_weights = {
            W_DELAY_init,
            W_CONFL_init,
            W_STRETCH_init,
            W_REC_init
    };

    mutable std::mt19937 metric_rng;
    int selectMetricIndex() const;
    double adaptive_heuristics_reaction_factor = 1;
    double adaptive_heuristics_decay_factor = 0.05;

    pair<int,int> roundRobin();
    pair<int, int> findMostDelayedAgentAndTime();

    enum SatHeuristic {
        SAT_ROUND_ROBIN = 0,
        SAT_MOST_DELAYED,
        SAT_ADAPTIVE
    };

    SatHeuristic sat_heuristic = SAT_ADAPTIVE;
    int          sat_submap_side = 5;
    int sat_prob_percent = 100;
    
    /* aggregated wall-clock times (seconds) */
    double sat_runtime_total   = 0.0;   ///< time spent in runSAT + immediate conflict-repair
    double other_runtime_total = 0.0;   ///< time spent in non-SAT operators

    destroy_heuristic fallback_destroy_strategy = INTERSECTION;   // used when SAT not selected
    std::string       fallback_replan_algo     = "PP";            // PP / CBS / EECBS

    /* helper function for converting string to enum from argument*/
    static destroy_heuristic strToDestroyHeuristic(const std::string& s)
    {
        if (s == "Random" || s == "RandomAgents")  return RANDOMAGENTS;
        if (s == "RandomWalk")                     return RANDOMWALK;
        if (s == "Intersection")                   return INTERSECTION;
        if (s == "SAT")                            return SAT;
        if (s == "Adaptive" || s == "ALNS")        return NONE;

        // default
        std::cerr << "[WARN] Unknown destroy heuristic '" << s
                  << "', falling back to Intersection.\n";
        return INTERSECTION;
    }
};
