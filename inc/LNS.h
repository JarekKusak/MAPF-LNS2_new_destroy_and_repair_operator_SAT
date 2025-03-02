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
    pair<vector<vector<int>>, vector<int>> getSubmapAndAgents(int agent_id, int submap_size, int agent_location); // helper function for getting sub-map
    vector<int> getAgentsToReplan(const vector<int>& agents_in_submap,
                                       const unordered_set<int>& submap_set,
                                       int problematic_timestep);
    void initializeSubmapData(const vector<vector<int>>& submap,
                                   unordered_set<int>& submap_set,
                                   unordered_map<int, pair<int, int>>& global_to_local);
    vector<vector<int>> generateMapRepresentation(const vector<vector<int>>& submap,
                                                       const vector<int>& agents_in_submap,
                                                       int problematic_timestep);
    bool solveWithSAT(vector<vector<int>>& map,
                           vector<pair<int, int>>& start_positions,
                           vector<pair<int, int>>& goal_positions,
                           vector<int>& agents_to_replan,
                           const vector<vector<int>>& submap,
                           int T_sync);
    void findStartAndGoalPositions(const vector<int>& agents_to_replan,
                                        const unordered_set<int>& submap_set,
                                        const unordered_map<int, pair<int, int>>& global_to_local,
                                        vector<pair<int, int>>& start_positions,
                                        vector<pair<int, int>>& goal_positions,
                                        int T_sync);
    int findSyncTimeAndEntryTimes(const vector<int>& agents_to_replan,
                                       const unordered_set<int>& submap_set,
                                       unordered_map<int, int>& agent_entry_time);
    void synchronizeAgentPaths(vector<int>& agents_to_replan,
                                    unordered_map<int, int>& agent_entry_time,
                                    int T_sync);
};
