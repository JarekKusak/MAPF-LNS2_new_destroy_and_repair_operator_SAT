#pragma once
#include "common.h"
#include "SpaceTimeAStar.h"
#include "SIPP.h"

#include <unordered_set>
#include <unordered_map>

class PathTable; // forward declaration

struct Agent
{
    int id;
    SingleAgentSolver* path_planner = nullptr; // start, goal, and heuristics are stored in the path planner
    Path path;

    Agent(const Instance& instance, int id, bool sipp) : id(id)
    {
        if(sipp)
            path_planner = new SIPP(instance, id);
        else
            path_planner = new SpaceTimeAStar(instance, id);
    }
    ~Agent(){ delete path_planner; }

    int getNumOfDelays() const
    {
        return (int) path.size() - 1 - path_planner->my_heuristic[path_planner->start_location];
    }

    /* analysis path of the agent and returns maximum delay count and the most problematic time step */
    pair<int,int> getMostProblematicDelay(
            const PathTable& path_table,
            const std::set<std::pair<int,int>>& ignored_ts) const {
        int max_delays = 0;
        int problematic_timestep = -1;

        for (int t = 0; t < path.size(); t++) {
            if (ignored_ts.count({id, t}))  // skip ignored
                continue;

            int delays = path_planner->getNumOfDelaysAtTimestep(
                    path_table, path, path[t].location, t);

            if (delays > max_delays) {
                max_delays = delays;
                problematic_timestep = t;
            }
        }
        return {max_delays, problematic_timestep};
    }

    vector<int> getTopDelayedTimesteps(const PathTable& tab,
                                       const set<pair<int,int>>& ignored) const {
        int maxDel = 0;
        vector<int> best;

        for (int t = 0; t < path.size(); ++t) {
            if (ignored.count({id,t})) continue;
            int d = path_planner->getNumOfDelaysAtTimestep(tab, path,
                                                           path[t].location, t);
            if (d == 0) continue;
            if (d > maxDel) { // new max
                maxDel = d;
                best.clear();
                best.push_back(t);
            } else if (d == maxDel) {
                best.push_back(t);
            }
        }
        return best; // empty = agent without delays...
    }

    struct AgentStats {
        int prev_delay_max = 0;
        int prev_conflict_cnt = 0;
        double prev_stretch_ratio = 0;
        int prev_last_replanned = 0;

        int   delay_max        = 0;
        int   conflict_cnt     = 0;
        double stretch_ratio   = 0.0;
        int   last_replanned   = 0;
    };
    mutable AgentStats stats;
};

struct Neighbor
{
    vector<int> agents;
    int sum_of_costs;
    int old_sum_of_costs;
    set<pair<int, int>> colliding_pairs;  // id1 < id2
    set<pair<int, int>> old_colliding_pairs;  // id1 < id2
    vector<Path> old_paths;

    // NOVÉ atributy pro SAT
    std::vector<vector<int>> submap;
    std::unordered_set<int> submap_set;
    std::unordered_map<int, pair<int, int>> global_to_local;
    int T_sync = -1;
    int key_agent_id = -1;
    std::vector<vector<int>> map;
};

class BasicLNS
{
public:
    // statistics
    int num_of_failures = 0; // #replanning that fails to find any solutions
    list<IterationStats> iteration_stats; // stats about each iteration
    double runtime = 0;
    double average_group_size = -1;
    int sum_of_costs = 0;

    BasicLNS(const Instance& instance, double time_limit, int neighbor_size, int screen);
    virtual string getSolverName() const = 0;
protected:
    // input params
    const Instance& instance; // avoid making copies of this variable as much as possible
    double time_limit;
    double replan_time_limit; // time limit for replanning
    int neighbor_size;
    int screen;

    // adaptive LNS
    bool ALNS = false;
    double decay_factor = -1;
    double reaction_factor = -1;
    vector<double> destroy_weights;
    int selected_neighbor;

    // helper variables
    high_resolution_clock::time_point start_time;
    Neighbor neighbor;

    void rouletteWheel();
};