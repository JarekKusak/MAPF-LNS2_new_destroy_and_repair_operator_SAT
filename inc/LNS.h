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

enum destroy_heuristic { RANDOMAGENTS, RANDOMWALK, INTERSECTION, DESTORY_COUNT, SAT };

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
        int screen, PIBTPPS_option pipp_option);
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

    pair<int, int> findMostDelayedAgent(); // the function is now returning the most problematic agent and his most problematic time step
    int findRandomAgent() const;
    void randomWalk(int agent_id, int start_location, int start_timestep,
                    set<int>& neighbor, int neighbor_size, int upperbound);

    bool generateNeighborBySAT(); // new destroy operator
    bool runSAT(); // new repair operator
    pair<vector<vector<int>>, vector<int>> getSubmapAndAgents(int agent_id, int submap_size, int agent_location, int timestep); // helper function for getting sub-map
    int last_selected_agent = -1;   // udržuje, na kom jsme skončili
    int countConflicts(const Agent& ag);
    int computeMaxDelay(const Agent& ag);
    double agentScore(const Agent& ag, double w_delay,
                      double w_conf, double w_stretch,
                      double w_recency) const;
    pair<int,int> findBestAgentAndTime();
    const double W_DELAY   = 4.0;
    const double W_CONFL   = 2.0;
    const double W_STRETCH = 1.0;
    const double W_REC     = 0.5;
    void updateAllStats(int iter);
    int current_iter;

};
