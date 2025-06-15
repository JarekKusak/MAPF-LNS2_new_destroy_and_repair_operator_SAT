#include "LNS.h"
#include "ECBS.h"
#include <queue>
#include <memory>

#include "../include/MAPF.hpp"
#include "SATUtils.h"
#include "Log.h"

LNS::LNS(const Instance& instance, double time_limit, const string & init_algo_name, const string & replan_algo_name,
         const string & destory_name, int neighbor_size, int num_of_iterations, bool use_init_lns,
         const string & init_destory_name, bool use_sipp, int screen, PIBTPPS_option pipp_option,
         const string& sat_heur_name, int sat_submap_side, int sat_prob_percent, const std::string& fallback_dest_name,
         const std::string& fallback_replan_name) :
         BasicLNS(instance, time_limit, neighbor_size, screen),
         init_algo_name(init_algo_name),  replan_algo_name(replan_algo_name),
         num_of_iterations(num_of_iterations),
         use_init_lns(use_init_lns),init_destory_name(init_destory_name),
         path_table(instance.map_size), pipp_option(pipp_option),
         sat_submap_side(sat_submap_side), sat_prob_percent(sat_prob_percent), fallback_replan_algo(fallback_replan_name),
         fallback_destroy_strategy(strToDestroyHeuristic(fallback_dest_name)) {
    start_time = Time::now();
    replan_time_limit = time_limit / 100;
    component_weights = {W_DELAY_init, W_CONFL_init, W_STRETCH_init, W_REC_init};
    if (destory_name == "Adaptive")
    {
        ALNS = true;
        destroy_weights.assign(DESTORY_COUNT, 1);
        decay_factor = 0.01;
        reaction_factor = 0.01;
    }
    else if (destory_name == "RandomWalk")
        destroy_strategy = RANDOMWALK;
    else if (destory_name == "Intersection")
        destroy_strategy = INTERSECTION;
    else if (destory_name == "Random")
        destroy_strategy = RANDOMAGENTS;
    else if (destory_name == "SAT") { // new operator
        destroy_strategy = SAT;
        //adaptive_heuristics_reaction_factor = 100;
        //adaptive_heuristics_decay_factor = 0.01;
    }
    else
    {
        cerr << "Destroy heuristic " << destory_name << " does not exists. " << endl;
        exit(-1);
    }

    if      (sat_heur_name == "roundRobin")  sat_heuristic = SAT_ROUND_ROBIN;
    else if (sat_heur_name == "mostDelayed") sat_heuristic = SAT_MOST_DELAYED;
    else                                     sat_heuristic = SAT_ADAPTIVE;   // default

    //fallback_destroy_strategy = strToDestroyHeuristic(fallback_dest_name);

    if (fallback_destroy_strategy == NONE)
    {
        SAT_DBG("ALNS set as fallback destroy strategy.");
        ALNS = true;
        destroy_weights.assign(DESTORY_COUNT, 1);
        decay_factor     = 0.01;
        reaction_factor  = 0.01;
    }

    int N = instance.getDefaultNumberOfAgents();
    agents.reserve(N);
    for (int i = 0; i < N; i++)
        agents.emplace_back(instance, i, use_sipp);
    preprocessing_time = ((fsec)(Time::now() - start_time)).count();
    if (screen >= 2)
        cout << "Pre-processing time = " << preprocessing_time << " seconds." << endl;
}

/* getting the submap around one or more agents and identifying agents in these submaps */
pair<vector<vector<int>>, vector<int>> LNS::getSubmapAndAgents(int agent_id, int submap_size, int agent_location, int timestep) {
    int map_width = instance.num_of_cols;
    int map_height = instance.num_of_rows;

    int submap_side = sqrt(submap_size); // assuming submap_size is a perfect square
    if (submap_side * submap_side != submap_size) {
        cout << "Error: Submap size is not a perfect square!" << endl;
        return {{}, {}};
    }

    // 2D submap initialized to -1
    vector<vector<int>> submap(submap_side, vector<int>(submap_side, -1));
    vector<int> agents_in_submap;
    set<int> conflicting_agents;

    // center of the submap
    int agent_x = agent_location / map_width; // most delayed agent is the center of submap
    int agent_y = agent_location % map_width;

    int half_side = submap_side / 2;

    for (int dx = -half_side; dx <= half_side; ++dx) {
        for (int dy = -half_side; dy <= half_side; ++dy) {
            int x = agent_x + dx; // relative displacement from the agent along the horizontal
            int y = agent_y + dy; // relative displacement from the agent along the vertical

            // ensure we are within map boundaries (if not -> skip)
            if (x >= 0 && x < map_height && y >= 0 && y < map_width) {
                int global_pos = x * map_width + y; // the unique index of a cell in the global map

                // map (dx, dy) to submap indices
                int submap_x = dx + half_side;
                int submap_y = dy + half_side;

                if (submap_x >= 0 && submap_x < submap_side && submap_y >= 0 && submap_y < submap_side) {
                    submap[submap_x][submap_y] = global_pos;
                    path_table.get_agents_at_timestep(conflicting_agents, global_pos, timestep);
                    //path_table.get_agents(conflicting_agents, global_pos); // collect all agents in the submap (at EVERY timestep)
                }
            }
        }
    }

    agents_in_submap.assign(conflicting_agents.begin(), conflicting_agents.end());
    return {submap, agents_in_submap};
}

// --------------------------------------------------------
// DESTROY phase: generateNeighborBySAT() – finds submap, agents, T_sync, etc.
// --------------------------------------------------------
bool LNS::generateNeighborBySAT() {
    SAT_DBG("====================");
    SAT_DBG("SAT destroy operator called.");

    pair<int,int> agent_time;
    switch (sat_heuristic) {
        case SAT_ROUND_ROBIN:    agent_time = roundRobin();                  break;
        case SAT_MOST_DELAYED:   agent_time = findMostDelayedAgentAndTime(); break;
        case SAT_ADAPTIVE:       agent_time = findBestAgentAndTime();        break;
    }
    SAT_DBG("SAT heuristic: " << sat_heuristic);
    auto [key_agent_id, problematic_timestep] = agent_time;
    if (key_agent_id < 0) {
        SAT_DBG("No delayed agent found.");
        ignored_agents_with_timestep.clear();
        return false;
    }
    SAT_DBG("key_agent_id: " << key_agent_id);
    SAT_DBG("key_agent_id global path length: " << agents[key_agent_id].path.size());

    int agent_loc = agents[key_agent_id].path[problematic_timestep].location;
    int submap_size = sat_submap_side * sat_submap_side;

    auto [submap, agents_in_submap] = getSubmapAndAgents(key_agent_id, submap_size, agent_loc, problematic_timestep);

    std::unordered_set<int> submap_set;
    std::unordered_map<int, pair<int, int>> global_to_local;
    SATUtils::initializeSubmapData(submap, submap_set, global_to_local);

    std::vector<vector<int>> map = SATUtils::generateMapRepresentation(submap, agents_in_submap, problematic_timestep, instance, agents);

    std::vector<int> agents_to_replan = SATUtils::getAgentsToReplan(agents_in_submap, submap_set, problematic_timestep, agents);
    if (agents_to_replan.empty()) {
        std::cout << "No agents to replan in submap." << std::endl;
        return false;
    }

    int T_sync = problematic_timestep;
    // Debug output - can be useful for diagnostics
    vector<pair<int,int>> start_positions, goal_positions;
    for (int agent : agents_to_replan) {
        int start_global = -1, goal_global = -1;
        int goal_time = -1;
        if ((size_t)T_sync >= agents[agent].path.size()) {
            std::cout << "Agent " << agent << " has no position defined at time T_sync!" << std::endl;
            continue;
        }
        int loc_at_Tsync = agents[agent].path[T_sync].location;
        if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
            std::cout << "Agent " << agent << " is not in the submap at time T_sync!" << std::endl;
            continue;
        }
        start_global = loc_at_Tsync;
        for (int t = T_sync; t < (int)agents[agent].path.size(); t++) {
            int location = agents[agent].path[t].location;
            if (submap_set.find(location) != submap_set.end()) {
                goal_global = location;
                goal_time = t;
            } else break;
        }
        if (goal_time == -1) {
            std::cout << "Agent " << agent << " has no valid goal position!" << std::endl;
            continue;
        }
        auto itS = global_to_local.find(start_global);
        auto itG = global_to_local.find(goal_global);
        if (itS == global_to_local.end() || itG == global_to_local.end()) {
            std::cout << "Missing mapping from global to local coordinates!" << std::endl;
            continue;
        }
        start_positions.push_back(itS->second);
        goal_positions.push_back(itG->second);
        SAT_DBG("Agent " << agent << " | Start (global): " << start_global
                 << " -> (local): (" << itS->second.first << ", " << itS->second.second << ")"
                 << " at time " << T_sync << " | Goal (global): " << goal_global
                 << " -> (local): (" << itG->second.first << ", " << itG->second.second << ")"
                 << " at time " << goal_time);
    }

    neighbor.agents = agents_to_replan;
    neighbor.submap = submap;
    neighbor.submap_set = submap_set;
    neighbor.global_to_local = global_to_local;
    neighbor.map = map;
    neighbor.T_sync = T_sync;
    neighbor.key_agent_id = key_agent_id;
    ignored_agents_with_timestep.insert({key_agent_id, T_sync});

    return true;
}

int LNS::computeMaxDelay(const Agent& ag)
{
    int best = 0;
    for (int t = 0; t < (int)ag.path.size(); ++t) {
        int d = ag.path_planner->getNumOfDelaysAtTimestep(
                path_table, ag.path, ag.path[t].location, t);
        best = std::max(best, d);
    }
    return best;
}

int LNS::countConflicts(const Agent& ag) {
    int c = 0;
    for (int t = 1; t < ag.path.size(); ++t)
        if (!instance.validMove(ag.path[t-1].location, ag.path[t].location))
            ++c;
    return c;
}

// The higher the score, the more "interesting" the agent is for further destruction/replan.
double LNS::agentScore(const Agent& ag) const {
    const auto& st = ag.stats;
    return  component_weights[0] * st.delay_max +
            component_weights[1] * st.conflict_cnt +
            component_weights[2] * st.stretch_ratio +
            component_weights[3] * (current_iter - st.last_replanned);
}

/*
    1. Find the agent with the highest agentScore(ag).
    2. For that agent, take its delay_max (maximum delay) and find all times t where this delay actually occurred.
    3. From all these time steps, randomly choose one that you haven't ignored yet (this must be tracked via the ignored_agents_with_timestep set).
    4. If all are ignored, clear the set and repeat the search; then if still nothing, return {-1,-1} -> do nothing.
    The result is the pair (key_agent_id, problematic_timestep).
 */
pair<int,int> LNS::findBestAgentAndTime()
{
    int best_id = -1;
    double best_score = -1;

    for (const auto& ag : agents) {
        double s = agentScore(ag);
        if (s < 0) continue;  // (could be deleted)

        if (s > best_score) {
            best_score = s;
            best_id = ag.id;
        }
    }
    if (best_id == -1) return {-1,-1};

    /* ‑‑ select timestep from largest delays ‑‑ */
    vector<int> ts;
    int dmax = agents[best_id].stats.delay_max;

    /* We find all time steps t where the delay is equal to dmax and which are not yet in the set ignored_agents_with_timestep. */
    for (int t = 0; t < (int)agents[best_id].path.size(); ++t) {
        int d = agents[best_id].path_planner->getNumOfDelaysAtTimestep(
                path_table,
                agents[best_id].path,
                agents[best_id].path[t].location,
                t);
        if (d == dmax && !ignored_agents_with_timestep.count({best_id, t}))
            ts.push_back(t);
    }

    /*
      If the first pass did not generate any t, it means that we have already “ignored” all of its worst times once.
      We clear the entire ignored_agents_with_timestep (but only once!) and repeat the selection once.
      If we still do not find any t, we return {-1,-1}, i.e. no action.
     */
    if (ts.empty()) {
        // all candidate timesteps were ignored; clear set and retry once
        ignored_agents_with_timestep.clear();
        // recompute timesteps for best_id
        ts.clear();
        for (int t = 0; t < (int)agents[best_id].path.size(); ++t) {
            int d = agents[best_id].path_planner->getNumOfDelaysAtTimestep(
                    path_table,
                    agents[best_id].path,
                    agents[best_id].path[t].location,
                    t);
            if (d == dmax && !ignored_agents_with_timestep.count({best_id, t}))
                ts.push_back(t);
        }
    }
    /* I called it recursively before
    if (ts.empty()) { // all of dmax are in the ignored list
        ignored_agents_with_timestep.clear(); // => allow new sequential pass
        return findBestAgentAndTime();   // recursion
    }*/

    if (ts.empty()) {
        // still no valid timestep
        return {-1, -1};
    }

    int chosen_t = ts[rand()%ts.size()];

    last_selected_agent = best_id;
    return {best_id, chosen_t};
}

void LNS::updateAllStats(int iter) {
    for (auto& ag : agents) {
        auto& st = ag.stats;

        // save previous values
        st.prev_delay_max     = st.delay_max;
        st.prev_conflict_cnt  = st.conflict_cnt;
        st.prev_stretch_ratio = st.stretch_ratio;
        st.prev_last_replanned = st.last_replanned;

        // maximum delay at any step of the path
        st.delay_max     = computeMaxDelay(ag);
        // number of invalid steps (edge/vertex conflicts)
        st.conflict_cnt  = countConflicts(ag);
        // relative path extension compared to the heuristic lower bound
        st.stretch_ratio = double(ag.path.size() - 1) /
                           ag.path_planner->my_heuristic[ ag.path_planner->start_location ];
    }
}

/*
   We want the metric that led to success (here index 0 = delay_max) to grow, the others to decay weakly:
   reaction_factor (e.g. 0.01) determines how strongly a weight is "rewarded" for a given delta.
   decay_factor (e.g. 0.01) determines how quickly the other weights return to greater versatility when they are not selected.
*/
void LNS::updateComponentWeights(int metric_index, double delta) {
    // Update selected metric weight multiplicatively and decay others
    for (int i = 0; i < (int)component_weights.size(); ++i) {
        if (i == metric_index)
            component_weights[i] *= (1.0 + adaptive_heuristics_reaction_factor * delta);
        else component_weights[i] *= (1.0 - adaptive_heuristics_decay_factor);
    }
    double sum_w = std::accumulate(component_weights.begin(), component_weights.end(), 0.0);
    if (sum_w > 0)
        for (auto &w : component_weights)
            w /= sum_w;
}

// --------------------------------------------------------
// REPAIR phase: runSAT() – calls findLocalPaths + solveWithSAT,
//               and updates agent paths + path_table
// --------------------------------------------------------
bool LNS::runSAT()
{
    SAT_DBG("[REPAIR] SAT destroy operator called.");
    SAT_DBG("[REPAIR] SAT operator – launching subproblem FOR OPTIMIZATION.");

    const auto& agents_to_replan = neighbor.agents;
    const auto& submap           = neighbor.submap;
    const auto& submap_set       = neighbor.submap_set;
    const auto& global_to_local  = neighbor.global_to_local;
    auto& map              = neighbor.map;
    int T_sync = neighbor.T_sync;
    int key_agent_id = neighbor.key_agent_id;

    auto local_paths = SATUtils::findLocalPaths(agents_to_replan, submap, submap_set, global_to_local, T_sync, agents);

    bool success = SATUtils::solveWithSAT(map, local_paths, agents_to_replan, submap, T_sync, agents);

    if (!success) {
        for (auto a : agents_to_replan)
            path_table.insertPath(agents[a].id, agents[a].path); // return old paths of agents

        SAT_STAT("SAT solver failed to find a valid solution.");
        return false;
    }
    else { // if succesful replan, add current_iter to agents' stats (for adaptive heur)
        for (auto a: agents_to_replan)
            agents[a].stats.last_replanned = current_iter;
    }

    // Additional validity check for paths returned by the SAT solver:
    // If an agent performs an invalid move (e.g. diagonal), revert to the old path.
    for (size_t idx = 0; idx < agents_to_replan.size(); ++idx) {
        int ag = agents_to_replan[idx];
        bool invalid_move = false;

        for (size_t t = 1; t < agents[ag].path.size(); ++t) {
            int from = agents[ag].path[t - 1].location;
            int to   = agents[ag].path[t].location;
            if (!instance.validMove(from, to)) {
                invalid_move = true;
                std::cout << "Agent " << ag << " made an invalid move " << from << " -> " << to
                                  << " between times " << t - 1 << " and " << t << ". Reverting to previous path." << std::endl;
                break;
            }
        }

        if (invalid_move)
            agents[ag].path = neighbor.old_paths[idx];
    }

    neighbor.sum_of_costs = 0;
    for (int ag : agents_to_replan)
        neighbor.sum_of_costs += (int)agents[ag].path.size() - 1;

    if (neighbor.sum_of_costs <= neighbor.old_sum_of_costs) {
        for (int a : agents_to_replan) {
            path_table.insertPath(agents[a].id, agents[a].path);
            SAT_DBG("(LNS.cpp) New path for agent " << a << ": ");
            {
                std::stringstream ss;
                for (auto loc: agents[a].path)
                    ss << loc.location << ", ";
                SAT_DBG(ss.str());
            }
        }

        auto& key_stats = agents[key_agent_id].stats;
        key_stats.prev_delay_max     = key_stats.delay_max;
        key_stats.prev_conflict_cnt  = key_stats.conflict_cnt;
        key_stats.prev_stretch_ratio = key_stats.stretch_ratio;

        key_stats.delay_max     = computeMaxDelay(agents[key_agent_id]);
        key_stats.conflict_cnt  = countConflicts(agents[key_agent_id]);
        key_stats.stretch_ratio = double(agents[key_agent_id].path.size()-1) /
                                  agents[key_agent_id].path_planner
                                          ->my_heuristic[agents[key_agent_id]
                                          .path_planner->start_location];

        double delta = double(neighbor.old_sum_of_costs - neighbor.sum_of_costs)
                       / double(neighbor.old_sum_of_costs);

        // how did statistics change for key agent
        const Agent& key_ag = agents[neighbor.key_agent_id];
        const auto& st = key_ag.stats;

        double d_delta = double(st.prev_delay_max    - st.delay_max); // drop  ==> positive
        double c_delta = double(st.prev_conflict_cnt - st.conflict_cnt);
        double s_delta =         st.prev_stretch_ratio - st.stretch_ratio; // drop  ==> positive
        // this ensures that the recency component gets a positive number only when the agent moves forward
        double l_delta = (delta>0) ? (st.last_replanned - st.prev_last_replanned) : 0; // only if SoC dropped

        // find the index of the component with the largest *relative* contribution (>0)
        std::array<double,4> deltas = {d_delta, c_delta, s_delta, l_delta};
        int metric_index = std::distance(deltas.begin(),
                                         std::max_element(deltas.begin(), deltas.end()));
        if (deltas[metric_index] <= 0)
            metric_index = 0; // fallback – reward delay when nothing else has improved

        // reward chosen weight
        updateComponentWeights(metric_index, delta);

        SAT_DBG("component_weights = {"
                        << component_weights[0] << ", "
                        << component_weights[1] << ", "
                        << component_weights[2] << ", "
                        << component_weights[3] << "}");

        return true;
    } else {
        SAT_DBG("[INFO] New SAT solution is worse, reverting.");
        for (int i = 0; i < (int)neighbor.agents.size(); i++) {
            int a = neighbor.agents[i];
            agents[a].path = neighbor.old_paths[i];
            SAT_DBG("(LNS.cpp) Reverted path for agent " << a << ": ");
            for (auto loc : agents[a].path)
                SAT_DBG(loc.location << ", ");

            path_table.insertPath(agents[a].id, agents[a].path);
            SAT_DBG("Verifying OLD path_table for agent " << a << ":");
            for (int t = 0; t < (int) agents[a].path.size(); t++)
            {
                int loc = agents[a].path[t].location;
                if (loc >= 0 && loc < (int) path_table.table.size())
                {
                    if ((int) path_table.table[loc].size() > t)
                        SAT_DBG("  time=" << t << ", loc=" << loc
                                 << ", table=" << path_table.table[loc][t]);
                    else
                        SAT_DBG("  time=" << t << ", loc=" << loc << " => out of range");
                }
            }
        }
        neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        return false;
    }
}

void LNS::rollbackNeighbor()
{
    if (!iter_backup_valid) return;        // nothing to return

    // return snapshot
    path_table.reset();

    for (int id = 0; id < agents.size(); ++id) {
        agents[id].path = iter_backup_paths[id];
        path_table.insertPath(id, agents[id].path);
    }

    sum_of_costs      = iter_backup_soc;
    iter_backup_valid = false; // snapshot used

    SAT_DBG("Rollback -> SoC = " << sum_of_costs << " (restored entire snapshot)");
}

// ============================================
// Helper function for immediate conflict repair via init_lns
// ============================================
void LNS::doInitLNSRepair(const string& debug_reason) {
    auto repair_start = Time::now();

    SAT_DBG("Attempting immediate repair via init_lns " << debug_reason << ".");
    // if there is no time left, we add 100 ms for the possibility of corrections so that the program does not fail validation
    //double repl_budget = std::max(0.1, time_limit - runtime);
    init_lns = new InitLNS(instance, agents, time_limit - runtime,
                           "PP",//replan_algo_name,
                           init_destory_name,
                           neighbor_size, screen);

    SAT_DBG("Passing " << agents.size()
                       << " agents to init_lns (skip=true).");

    // ------------------------------------------------------------------
    // Output: path_table contents for selected agents (e.g., neighbor.agents)
    // ------------------------------------------------------------------
    SAT_DBG("path_table contents for selected agents (neighbor.agents):");
    for (int a : neighbor.agents) {
        SAT_DBG("  Agent " << a << " => controlling path length=" << agents[a].path.size());
        for (int t = 0; t < (int)agents[a].path.size(); t++) {
            int loc = agents[a].path[t].location;
            // Check if loc is a valid index
            if (loc < 0 || loc >= (int)path_table.table.size()) {
                SAT_DBG("    [time=" << t << "]: loc=" << loc << " (out of range)");
                continue;
            }
            // Check if path_table.table[loc].size() > t
            if ((int)path_table.table[loc].size() <= t) {
                SAT_DBG("    [time=" << t << ", loc=" << loc << "]: path_table.table[loc].size()="
                                     << path_table.table[loc].size() << " => out of range for t=" << t);
                continue;
            }
        }
    }

    bool fixed = init_lns->run(true);

    SAT_DBG("init_lns->sum_of_costs after init_lns->run: " << init_lns->sum_of_costs);
    if (fixed) {
        sum_of_costs = init_lns->sum_of_costs;
        //neighbor.old_sum_of_costs = init_lns->sum_of_costs;

        SAT_DBG("sum_of_costs after assignment from init_lns->run: " << sum_of_costs);
        for (const auto &agent : agents)
            path_table.insertPath(agent.id, agent.path);

        init_lns->clear();
    }
    else {
        std::cout << "[ERROR] Could not repair solution right after SAT." << std::endl;
        rollbackNeighbor();
    }
    other_runtime_total += ((fsec)(Time::now() - repair_start)).count();
}

bool LNS::run()
{
    // Open file for logging output
    std::ofstream out("log.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(out.rdbuf());

    sat_runtime_total   = 0.0;
    other_runtime_total = 0.0;
    overhead_runtime_total = 0.0;

    sum_of_distances = 0;
    for (const auto & agent : agents)
        sum_of_distances += agent.path_planner->my_heuristic[agent.path_planner->start_location];

    initial_solution_runtime = 0;
    start_time = Time::now();
    bool succ = getInitialSolution();
    initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();

    if (!succ && initial_solution_runtime < time_limit)
    {
        if (use_init_lns)
        {
            init_lns = new InitLNS(instance, agents, time_limit - initial_solution_runtime,
                                   "PP",//replan_algo_name,
                                   init_destory_name, neighbor_size, screen);
            succ = init_lns->run(false);
            if (succ)
            {
                path_table.reset();
                for (const auto & agent : agents)
                    path_table.insertPath(agent.id, agent.path);
                init_lns->clear();
                initial_sum_of_costs = init_lns->sum_of_costs;
                sum_of_costs = initial_sum_of_costs;
            }
            initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
        }
        else
        {
            while (!succ && initial_solution_runtime < time_limit)
            {
                succ = getInitialSolution();
                initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
                restart_times++;
            }
        }
    }

    iteration_stats.emplace_back(neighbor.agents.size(),
                                 initial_sum_of_costs, initial_solution_runtime, init_algo_name);

    other_runtime_total += initial_solution_runtime;
    runtime = initial_solution_runtime;
    if (!succ) {
        std::cout << "[ERROR] Failed to find an initial solution in "
                  << runtime << " seconds after  " << restart_times << " restarts" << std::endl;
        return false;
    }

    bool needConflictRepair = false;

    size_t sat_iter_count        = 0;     // how many times the SAT actually ran
    double last_sat_iter_runtime = 0.0;   // length of the last SAT iteration
    double last_other_iter_runtime = 0.0;
    double max_sat_iter_runtime  = 0.0;   // longest SAT iteration

    // Maximum SAT attempts per outer iteration
    const int MAX_SAT_TRIALS = 8;   // maximum SAT attempts within one outer iteration

    /* ------------------------------------------------------
    Flags fixing the choice of operator for the current
    outer iteration (i.e. before we successfully complete it
    or "continue" returns control to the beginning of the loop).
     ------------------------------------------------------ */
    bool decision_taken = false;   // Have we already drawn?
    bool SATchosenIter  = false;   // was SAT chosen for this iteration?

    // Optimization loop
    while (runtime < time_limit && iteration_stats.size() <= num_of_iterations) {
        current_iter = static_cast<int>(iteration_stats.size());
        auto iter_begin_TS = Time::now(); // framework overhead runtime
        selected_neighbor = -1;
        last_other_iter_runtime = 0.0;

        // one-time SAT x fallback selection for this iteration
        if (!decision_taken) {
            if (destroy_strategy == SAT)
                SATchosenIter = (rand() % 100) < sat_prob_percent;
            else SATchosenIter = false; // pure non-SAT configuration
            decision_taken = true;
        }

        SAT_STAT("Iteration " << current_iter);
        std::cout.flush();
        runtime = ((fsec)(Time::now() - start_time)).count();

        if (ALNS && destroy_strategy != SAT)
            chooseDestroyHeuristicbyALNS();

        bool opSuccess = false;
        bool SATchosen = false;

        updateAllStats(current_iter);

        if (destroy_strategy == SAT) { // SAT is destroy operator merged with replan operator
            if (SATchosenIter) {
                SATchosen = true;
                SAT_DBG("Using SAT operator (destroy+repair SAT).");

                auto sat_start = Time::now();

                // backup snapshot
                iter_backup_soc   = sum_of_costs;
                iter_backup_valid = true;
                iter_backup_paths.resize(agents.size());
                for (int i = 0; i < agents.size(); ++i)
                    iter_backup_paths[i] = agents[i].path; // full deep-copy - INEFFECTIVE - lot of consumption memory and time

                // New SAT trial loop with global time check
                for (int sat_trials = 0; !opSuccess && sat_trials < MAX_SAT_TRIALS &&
                 ((fsec)(Time::now() - start_time)).count() < time_limit;
                    ++sat_trials) {
                    if (!generateNeighborBySAT())
                        continue;

                    /* --------- backup old paths and empty path_table --------- */
                    neighbor.old_paths.resize(neighbor.agents.size());
                    neighbor.old_sum_of_costs = 0;
                    for (int i = 0; i < (int)neighbor.agents.size(); ++i) {
                        int a = neighbor.agents[i];
                        neighbor.old_paths[i] = agents[a].path;
                        path_table.deletePath(a, agents[a].path);
                        neighbor.old_sum_of_costs += (int)agents[a].path.size() - 1;
                    }

                    opSuccess = runSAT();
                }

                double sat_iter_runtime =
                        ((fsec)(Time::now() - sat_start)).count();

                sat_runtime_total += sat_iter_runtime;
                last_sat_iter_runtime = sat_iter_runtime;
                if (sat_iter_runtime > max_sat_iter_runtime)
                    max_sat_iter_runtime = sat_iter_runtime;
                ++sat_iter_count;
            }
            else SAT_STAT("Random chance did not select SAT operator, using default destroy strategy " << fallback_destroy_strategy << " with replan algo " << fallback_replan_algo);
            SAT_DBG("opSuccess value: " << opSuccess);
        }

        // Fallback and SAT-failed guards
        if (!opSuccess && !SATchosen)
        {
            auto other_start = Time::now();
            // fallback neighbor generation
            int DEFAULT_DESTROY_STRATEGY;
            std::string DEFAULT_REPLAN_ALGO;

            if (destroy_strategy == SAT && !SATchosenIter) { // mix SAT with other operator from argument (satProb 1-99)
                DEFAULT_DESTROY_STRATEGY = fallback_destroy_strategy;
                DEFAULT_REPLAN_ALGO      = fallback_replan_algo;
            }
            else { // PURE non-SAT
                DEFAULT_DESTROY_STRATEGY = destroy_strategy;  // command line arguments
                DEFAULT_REPLAN_ALGO = replan_algo_name;  // PP / CBS / EECBS ...
            }

            // --- fallback neighbour generation ---------------------------------
            auto saved_destroy_strategy = destroy_strategy;   // SAT
            if (DEFAULT_DESTROY_STRATEGY == NONE)
            {
                /*  Adaptive fallback:
                    – let ALNS choose the destruction heuristic,
                    – save it in DEFAULT_DESTROY_STRATEGY,
                    – and then everything works as it should in the code.
                 */
                chooseDestroyHeuristicbyALNS();
                DEFAULT_DESTROY_STRATEGY = destroy_strategy;   // value from ALNS
                destroy_strategy = saved_destroy_strategy;
                SAT_DBG("ALNS picked destroy = " << DEFAULT_DESTROY_STRATEGY);
            }

            SAT_DBG("running destroy strategy " << DEFAULT_DESTROY_STRATEGY);
            SAT_DBG("running replan strategy " << DEFAULT_REPLAN_ALGO);

            switch (DEFAULT_DESTROY_STRATEGY)
            {
                case RANDOMWALK:
                    opSuccess = generateNeighborByRandomWalk();
                    break;
                case INTERSECTION:
                    opSuccess = generateNeighborByIntersection();
                    break;
                case RANDOMAGENTS:
                {
                    neighbor.agents.resize(agents.size());
                    for (int i = 0; i < (int)agents.size(); i++) neighbor.agents[i] = i;
                    if ((int)neighbor.agents.size() > neighbor_size)
                    {
                        std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
                        neighbor.agents.resize(neighbor_size);
                    }
                    opSuccess = true;
                }
                    break;
                default:
                    std::cerr << "Wrong neighbor generation strategy" << std::endl;
                    exit(-1);
            }

            if (!opSuccess) {
                double iter_total     = ((fsec)(Time::now() - iter_begin_TS)).count();
                double iter_accounted = SATchosen ? last_sat_iter_runtime : 0.0;   // re-plan neběžel
                overhead_runtime_total += std::max(0.0, iter_total - iter_accounted);
                continue;
            }


            neighbor.old_paths.resize(neighbor.agents.size());
            neighbor.old_sum_of_costs = 0;
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                int a = neighbor.agents[i];
                neighbor.old_paths[i] = agents[a].path;
                path_table.deletePath(a, agents[a].path);
                neighbor.old_sum_of_costs += (int)agents[a].path.size() - 1;
            }

            if      (DEFAULT_REPLAN_ALGO == "PP")   succ = runPP();
            else if (DEFAULT_REPLAN_ALGO == "CBS")  succ = runCBS();
            else if (DEFAULT_REPLAN_ALGO == "EECBS") succ = runEECBS();
            else { std::cerr << "Wrong replanning strategy " << std::endl; exit(-1); }

            auto other_iter_runtime = ((fsec)(Time::now() - other_start)).count();
            other_runtime_total += other_iter_runtime;
            last_other_iter_runtime = other_iter_runtime;
        }
        else if (!opSuccess && SATchosen)
        {
            // SAT was the chosen operator but didn’t find a solution within the limits:
            // account overhead and continue with next outer iteration
            double iter_total = ((fsec)(Time::now() - iter_begin_TS)).count();
            overhead_runtime_total += iter_total;
            continue;
        }
        else succ = opSuccess;// opSuccess = true => runSAT completed

        if (!succ)
            continue;

        // ALNS evaluation
        if (ALNS && selected_neighbor >= 0)
        {
            if (neighbor.old_sum_of_costs > neighbor.sum_of_costs)
                destroy_weights[selected_neighbor] =
                        reaction_factor * (neighbor.old_sum_of_costs - neighbor.sum_of_costs) / neighbor.agents.size()
                        + (1 - reaction_factor) * destroy_weights[selected_neighbor];
            else
                destroy_weights[selected_neighbor] =
                        (1 - decay_factor) * destroy_weights[selected_neighbor];
        }

        double iter_total     = ((fsec)(Time::now() - iter_begin_TS)).count();
        double iter_accounted = SATchosen ? last_sat_iter_runtime : last_other_iter_runtime;
        overhead_runtime_total += std::max(0.0, iter_total - iter_accounted);

        runtime = ((fsec)(Time::now() - start_time)).count();

        SAT_STAT("neighbor.sum_of_costs before recomputation: " << neighbor.sum_of_costs);
        SAT_STAT("neighbor.old_sum_of_costs before recomputation: " << neighbor.old_sum_of_costs);
        SAT_STAT("sum_of_costs before recomputation: " << sum_of_costs);

        // ------------------------------------------------
        // After SAT => validation and possible conflict repair
        // ------------------------------------------------
        if (destroy_strategy == SAT && opSuccess && SATchosen)
        {
            /* ---  synchronize global sum_of_costs *before* validation  ---
               runSAT has just changed the paths of selected agents and filled
               neighbor.{old_,}sum_of_costs.  validateSolution() checks
               the consistency of sum_of_costs, so we must
               recalculate it here, otherwise it will lag by half a step and report
               “sum of costs mismatch”.                                          */
            sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;
            /* so that delta is not added again a few lines below */
            neighbor.old_sum_of_costs = neighbor.sum_of_costs;

            SAT_DBG("Validate solution immediately after SAT success.");
            try {
                validateSolution();
                SAT_DBG("No problems after SAT replan.");
            } catch (const ValidationException& e) {
                std::cout << "[WARNING] Problem after SAT: " << e.what() << std::endl;
                // unify
                doInitLNSRepair("because problem occurred after SAT (should be applied only for conflicts...)");
                validateSolution();
            }
        } else sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;

        SAT_DBG("Recomputing sum_of_cost by dividing neighbor.sum_of_costs and neighbor.old_sum_of_costs");
        SAT_STAT("sum_of_costs after recomputation: " << sum_of_costs);

        if (screen >= 1)
        {
            SAT_STAT("Iteration " << iteration_stats.size()
                 << ", group size = " << neighbor.agents.size()
                 << ", solution cost = " << sum_of_costs
                 << ", remaining time = " << time_limit - runtime
                 );
        }
        iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, replan_algo_name);

        // prepare another outer iteration – we will draw again
        decision_taken = false;
        SATchosenIter  = false;
    }

    average_group_size = -iteration_stats.front().num_of_agents;
    for (const auto& data : iteration_stats)
        average_group_size += data.num_of_agents;
    if (average_group_size > 0)
        average_group_size /= (double)(iteration_stats.size() - 1);

    SAT_STAT(getSolverName()
         << ": runtime = " << runtime
         << ", iterations = " << iteration_stats.size()
         << ", solution cost = " << sum_of_costs
         << ", initial solution cost = " << initial_sum_of_costs
         << ", failed iterations = " << num_of_failures
         );

    SAT_STAT("SAT total runtime       = " << sat_runtime_total      << " s");
    SAT_STAT("Other operators runtime = " << other_runtime_total    << " s");
    SAT_STAT("Framework overhead      = " << overhead_runtime_total << " s");

    double recog = sat_runtime_total + other_runtime_total + overhead_runtime_total;
    SAT_STAT("Sanity check - sum      = " << recog
                                          << " s  (whole run: " << runtime << " s)");
    if (runtime > 0)
        SAT_STAT("SAT runtime ratio    = "
                         << 100.0 * sat_runtime_total / runtime << " %");

    std::cout.rdbuf(coutbuf);
    return true;
}

bool LNS::getInitialSolution()
{
    neighbor.agents.resize(agents.size());
    for (int i = 0; i < (int)agents.size(); i++)
        neighbor.agents[i] = i;
    neighbor.old_sum_of_costs = MAX_COST;
    neighbor.sum_of_costs = 0;
    bool succ = false;
    if (init_algo_name == "EECBS")
        succ = runEECBS();
    else if (init_algo_name == "PP")
        succ = runPP();
    else if (init_algo_name == "PIBT")
        succ = runPIBT();
    else if (init_algo_name == "PPS")
        succ = runPPS();
    else if (init_algo_name == "winPIBT")
        succ = runWinPIBT();
    else if (init_algo_name == "CBS")
        succ = runCBS();
    else
    {
        cerr <<  "Initial MAPF solver " << init_algo_name << " does not exist!" << endl;
        exit(-1);
    }
    if (succ)
    {
        initial_sum_of_costs = neighbor.sum_of_costs;
        sum_of_costs = neighbor.sum_of_costs;
        return true;
    }
    else
    {
        return false;
    }

}

bool LNS::runEECBS()
{
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(agents[i].path_planner);
    }

    ECBS ecbs(search_engines, screen - 1, &path_table);
    ecbs.setPrioritizeConflicts(true);
    ecbs.setDisjointSplitting(false);
    ecbs.setBypass(true);
    ecbs.setRectangleReasoning(true);
    ecbs.setCorridorReasoning(true);
    ecbs.setHeuristicType(heuristics_type::WDG, heuristics_type::GLOBAL);
    ecbs.setTargetReasoning(true);
    ecbs.setMutexReasoning(false);
    ecbs.setConflictSelectionRule(conflict_selection::EARLIEST);
    ecbs.setNodeSelectionRule(node_selection::NODE_CONFLICTPAIRS);
    ecbs.setSavingStats(false);
    double w;
    if (iteration_stats.empty())
        w = 5; // initial run
    else
        w = 1.1; // replan
    ecbs.setHighLevelSolver(high_level_solver_type::EES, w);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime;
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = ecbs.solve(T, 0);
    if (succ && ecbs.solution_cost < neighbor.old_sum_of_costs) // accept new paths
    {
        auto id = neighbor.agents.begin();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *ecbs.paths[i];
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = ecbs.solution_cost;
        if (sum_of_costs_lowerbound < 0)
            sum_of_costs_lowerbound = ecbs.getLowerBound();
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id, agents[id].path);
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
        if (!succ)
            num_of_failures++;
    }
    return succ;
}
bool LNS::runCBS()
{
    if (screen >= 2)
        cout << "old sum of costs = " << neighbor.old_sum_of_costs << endl;
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(agents[i].path_planner);
    }

    CBS cbs(search_engines, screen - 1, &path_table);
    cbs.setPrioritizeConflicts(true);
    cbs.setDisjointSplitting(false);
    cbs.setBypass(true);
    cbs.setRectangleReasoning(true);
    cbs.setCorridorReasoning(true);
    cbs.setHeuristicType(heuristics_type::WDG, heuristics_type::ZERO);
    cbs.setTargetReasoning(true);
    cbs.setMutexReasoning(false);
    cbs.setConflictSelectionRule(conflict_selection::EARLIEST);
    cbs.setNodeSelectionRule(node_selection::NODE_CONFLICTPAIRS);
    cbs.setSavingStats(false);
    cbs.setHighLevelSolver(high_level_solver_type::ASTAR, 1);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime; // time limit
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = cbs.solve(T, 0);
    if (succ && cbs.solution_cost <= neighbor.old_sum_of_costs) // accept new paths
    {
        auto id = neighbor.agents.begin();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *cbs.paths[i];
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = cbs.solution_cost;
        if (sum_of_costs_lowerbound < 0)
            sum_of_costs_lowerbound = cbs.getLowerBound();
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id, agents[id].path);
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;

        }
        if (!succ)
            num_of_failures++;
    }
    return succ;
}
bool LNS::runPP()
{
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());
    if (screen >= 2) {
        for (auto id : shuffled_agents)
            cout << id << "(" << agents[id].path_planner->my_heuristic[agents[id].path_planner->start_location] <<
                "->" << agents[id].path.size() - 1 << "), ";
        cout << endl;
    }
    int remaining_agents = (int)shuffled_agents.size();
    auto p = shuffled_agents.begin();
    neighbor.sum_of_costs = 0;
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime; // time limit
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    auto time = Time::now();
    ConstraintTable constraint_table(instance.num_of_cols, instance.map_size, &path_table);
    while (p != shuffled_agents.end() && ((fsec)(Time::now() - time)).count() < T)
    {
        int id = *p;
        if (screen >= 3)
            cout << "Remaining agents = " << remaining_agents <<
                 ", remaining time = " << T - ((fsec)(Time::now() - time)).count() << " seconds. " << endl
                 << "Agent " << agents[id].id << endl;
        agents[id].path = agents[id].path_planner->findPath(constraint_table);
        if (agents[id].path.empty()) break;
        neighbor.sum_of_costs += (int)agents[id].path.size() - 1;
        if (neighbor.sum_of_costs >= neighbor.old_sum_of_costs)
            break;
        remaining_agents--;
        path_table.insertPath(agents[id].id, agents[id].path);
        ++p;
    }
    if (remaining_agents == 0 && neighbor.sum_of_costs <= neighbor.old_sum_of_costs) // accept new paths
    {
        return true;
    }
    else // stick to old paths
    {
        if (p != shuffled_agents.end())
            num_of_failures++;
        auto p2 = shuffled_agents.begin();
        while (p2 != p)
        {
            int a = *p2;
            path_table.deletePath(agents[a].id, agents[a].path);
            ++p2;
        }
        if (!neighbor.old_paths.empty())
        {
            p2 = neighbor.agents.begin();
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                int a = *p2;
                agents[a].path = neighbor.old_paths[i];
                path_table.insertPath(agents[a].id, agents[a].path);
                ++p2;
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
        return false;
    }
}
bool LNS::runPPS(){
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);
    P.setTimestepLimit(pipp_option.timestepLimit);

    // seed for solver
    auto* MT_S = new std::mt19937(0);
    PPS solver(&P,MT_S);
    solver.setTimeLimit(time_limit);
//    solver.WarshallFloyd();
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}
bool LNS::runPIBT(){
    auto shuffled_agents = neighbor.agents;
     std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);

    // seed for solver
    auto MT_S = new std::mt19937(0);
    PIBT solver(&P,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}
bool LNS::runWinPIBT(){
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);
    P.setTimestepLimit(pipp_option.timestepLimit);

    // seed for solver
    auto MT_S = new std::mt19937(0);
    winPIBT solver(&P,pipp_option.windowSize,pipp_option.winPIBTSoft,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}

MAPF LNS::preparePIBTProblem(vector<int>& shuffled_agents){

    // seed for problem and graph
    auto MT_PG = new std::mt19937(0);

//    Graph* G = new SimpleGrid(instance);
    Graph* G = new SimpleGrid(instance.getMapFile());

    std::vector<Task*> T;
    PIBT_Agents A;

    for (int i : shuffled_agents){
        assert(G->existNode(agents[i].path_planner->start_location));
        assert(G->existNode(agents[i].path_planner->goal_location));
        auto a = new PIBT_Agent(G->getNode( agents[i].path_planner->start_location));

//        PIBT_Agent* a = new PIBT_Agent(G->getNode( agents[i].path_planner.start_location));
        A.push_back(a);
        Task* tau = new Task(G->getNode( agents[i].path_planner->goal_location));


        T.push_back(tau);
        if(screen>=5){
            cout<<"Agent "<<i<<" start: " <<a->getNode()->getPos()<<" goal: "<<tau->getG().front()->getPos()<<endl;
        }
    }

    return MAPF(G, A, T, MT_PG);

}

void LNS::updatePIBTResult(const PIBT_Agents& A, vector<int>& shuffled_agents){
    int soc = 0;
    for (int i=0; i<A.size();i++){
        int a_id = shuffled_agents[i];

        agents[a_id].path.resize(A[i]->getHist().size());
        int last_goal_visit = 0;
        if(screen>=2)
            std::cout<<A[i]->logStr()<<std::endl;
        for (int n_index = 0; n_index < A[i]->getHist().size(); n_index++){
            auto n = A[i]->getHist()[n_index];
            agents[a_id].path[n_index] = PathEntry(n->v->getId());

            //record the last time agent reach the goal from a non-goal vertex.
            if(agents[a_id].path_planner->goal_location == n->v->getId()
                && n_index - 1>=0
                && agents[a_id].path_planner->goal_location !=  agents[a_id].path[n_index - 1].location)
                last_goal_visit = n_index;

        }
        //resize to last goal visit time
        agents[a_id].path.resize(last_goal_visit + 1);
        if(screen>=2)
            std::cout<<" Length: "<< agents[a_id].path.size() <<std::endl;
        if(screen>=5){
            cout <<"Agent "<<a_id<<":";
            for (auto loc : agents[a_id].path){
                cout <<loc.location<<",";
            }
            cout<<endl;
        }
        path_table.insertPath(agents[a_id].id, agents[a_id].path);
        soc += (int)agents[a_id].path.size()-1;
    }

    neighbor.sum_of_costs =soc;
}

void LNS::chooseDestroyHeuristicbyALNS()
{
    rouletteWheel();
    switch (selected_neighbor)
    {
        case 0 : destroy_strategy = RANDOMWALK; break;
        case 1 : destroy_strategy = INTERSECTION; break;
        case 2 : destroy_strategy = RANDOMAGENTS; break;
        default : cerr << "ERROR" << endl; exit(-1);
    }
}

bool LNS::generateNeighborByIntersection()
{
    if (intersections.empty())
    {
        for (int i = 0; i < instance.map_size; i++)
        {
            if (!instance.isObstacle(i) && instance.getDegree(i) > 2)
                intersections.push_back(i);
        }
    }

    set<int> neighbors_set;
    auto pt = intersections.begin();
    std::advance(pt, rand() % intersections.size());
    int location = *pt;
    path_table.get_agents(neighbors_set, neighbor_size, location);
    if (neighbors_set.size() < neighbor_size)
    {
        set<int> closed;
        closed.insert(location);
        std::queue<int> open;
        open.push(location);
        while (!open.empty() && (int) neighbors_set.size() < neighbor_size)
        {
            int curr = open.front();
            open.pop();
            for (auto next : instance.getNeighbors(curr))
            {
                if (closed.count(next) > 0)
                    continue;
                open.push(next);
                closed.insert(next);
                if (instance.getDegree(next) >= 3)
                {
                    path_table.get_agents(neighbors_set, neighbor_size, next);
                    if ((int) neighbors_set.size() == neighbor_size)
                        break;
                }
            }
        }
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (neighbor.agents.size() > neighbor_size)
    {
        std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
        neighbor.agents.resize(neighbor_size);
    }
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by intersection " << location << endl;
    return true;
}

bool LNS::generateNeighborByRandomWalk()
{
    if (neighbor_size >= (int)agents.size())
    {
        neighbor.agents.resize(agents.size());
        for (int i = 0; i < (int)agents.size(); i++)
            neighbor.agents[i] = i;
        return true;
    }

    int a = findMostDelayedAgent();

    if (a < 0)
        return false;
    
    set<int> neighbors_set;
    neighbors_set.insert(a);
    randomWalk(a, agents[a].path[0].location, 0, neighbors_set, neighbor_size, (int) agents[a].path.size() - 1);
    int count = 0;
    while (neighbors_set.size() < neighbor_size && count < 10)
    {
        int t = rand() % agents[a].path.size();
        randomWalk(a, agents[a].path[t].location, t, neighbors_set, neighbor_size, (int) agents[a].path.size() - 1);
        count++;
        // select the next agent randomly
        int idx = rand() % neighbors_set.size();
        int i = 0;
        for (auto n : neighbors_set)
        {
            if (i == idx)
            {
                a = n;
                break;
            }
            i++;
        }
    }
    if (neighbors_set.size() < 2)
        return false;
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by random walks of agent " << a
             << "(" << agents[a].path_planner->my_heuristic[agents[a].path_planner->start_location]
             << "->" << agents[a].path.size() - 1 << ")" << endl;

    return true;
}

int LNS::findMostDelayedAgent()
{
    int a = -1;
    int max_delays = -1;
    for (int i = 0; i < agents.size(); i++)
    {
        if (tabu_list.find(i) != tabu_list.end())
            continue;
        int delays = agents[i].getNumOfDelays();
        if (max_delays < delays)
        {
            a = i;
            max_delays = delays;
        }
    }
    if (max_delays == 0)
    {
        tabu_list.clear();
        return -1;
    }
    tabu_list.insert(a);
    if (tabu_list.size() == agents.size())
        tabu_list.clear();
    return a;
}


pair<int,int> LNS::roundRobin()
{
    const int N = static_cast<int>(agents.size());
    if (N == 0) return {-1,-1};

    // Projdeme agenty cyklicky: (last+1) … (last+N)
    for (int step = 1; step <= N; ++step)
    {
        int idx = (last_selected_agent + step) % N;
        const auto& ag = agents[idx];

        // najdi NEJvíc zpožděný timestep tohoto agenta,
        // který zatím není v ignored_agents_with_timestep
        auto [delays, ts] =
                ag.getMostProblematicDelay(path_table, ignored_agents_with_timestep);

        if (ts < 0) // už nemá nic nového
            continue;

        // máme kandidáta
        last_selected_agent = idx;          // zapamatuj si ho pro příště
        return {idx, ts};

        /*
        auto ts_vec = ag.getTopDelayedTimesteps(path_table, ignored_agents_with_timestep);
        if (ts_vec.empty()) continue;

        int ts = ts_vec[rand() % ts_vec.size()]; // náhodný z nejhorších
        last_selected_agent = idx;
        return {idx, ts};*/
    }

    // nikdo nemá další nepokrytý timestep
    return {-1,-1};
}

pair<int, int> LNS::findMostDelayedAgentAndTime() {
    int max_delays = -1;
    int agent_with_max_delays = -1;
    int most_problematic_timestep = -1;

    for (const auto& agent : agents) {
        // přeskoč agenty, kteří už selhali
        auto [agent_max_delays, problematic_timestep] =
                agent.getMostProblematicDelay(path_table, ignored_agents_with_timestep);

        if (ignored_agents_with_timestep.count({agent.id, problematic_timestep}))
            continue;    // tuto (agent,timestep) dvojici už jsme zkoušeli

        if (agent_max_delays > max_delays) {
            max_delays = agent_max_delays;
            agent_with_max_delays = agent.id;
            most_problematic_timestep = problematic_timestep;
        }
    }

    return {agent_with_max_delays, most_problematic_timestep};
}

int LNS::findRandomAgent() const
{
    int a = 0;
    int pt = rand() % (sum_of_costs - sum_of_distances) + 1;
    int sum = 0;
    for (; a < (int) agents.size(); a++)
    {
        sum += agents[a].getNumOfDelays();
        if (sum >= pt)
            break;
    }
    assert(sum >= pt);
    return a;
}

// a random walk with path that is shorter than upperbound and has conflicting with neighbor_size agents
void LNS::randomWalk(int agent_id, int start_location, int start_timestep,
                     set<int>& conflicting_agents, int neighbor_size, int upperbound)
{
    int loc = start_location;
    for (int t = start_timestep; t < upperbound; t++)
    {
        auto next_locs = instance.getNeighbors(loc);
        next_locs.push_back(loc);
        while (!next_locs.empty())
        {
            int step = rand() % next_locs.size();
            auto it = next_locs.begin();
            advance(it, step);
            int next_h_val = agents[agent_id].path_planner->my_heuristic[*it];
            if (t + 1 + next_h_val < upperbound) // move to this location
            {
                path_table.getConflictingAgents(agent_id, conflicting_agents, loc, *it, t + 1);
                loc = *it;
                break;
            }
            next_locs.erase(it);
        }
        if (next_locs.empty() || conflicting_agents.size() >= neighbor_size)
            break;
    }
}

void LNS::validateSolution() const
{
    int sum = 0;
    for (const auto& a1_ : agents)
    {
        if (a1_.path.empty())
        {
            throw ValidationException("No solution for agent " + std::to_string(a1_.id));
        }
        else if (a1_.path_planner->start_location != a1_.path.front().location)
        {
            throw ValidationException("The path of agent " + std::to_string(a1_.id) +
                                      " starts from location " + std::to_string(a1_.path.front().location) +
                                      ", which is different from its start location " +
                                      std::to_string(a1_.path_planner->start_location));
        }
        else if (a1_.path_planner->goal_location != a1_.path.back().location)
        {
            throw ValidationException("The path of agent " + std::to_string(a1_.id) +
                                      " ends at location " + std::to_string(a1_.path.back().location) +
                                      ", which is different from its goal location " +
                                      std::to_string(a1_.path_planner->goal_location));
        }

        for (int t = 1; t < (int)a1_.path.size(); t++)
        {
            if (!instance.validMove(a1_.path[t - 1].location, a1_.path[t].location))
            {
                throw ValidationException("The path of agent " + std::to_string(a1_.id) +
                                          " jumps from " + std::to_string(a1_.path[t - 1].location) +
                                          " to " + std::to_string(a1_.path[t].location) +
                                          " between timesteps " + std::to_string(t - 1) + " and " + std::to_string(t));
            }
        }

        sum += (int)a1_.path.size() - 1;

        for (const auto& a2_ : agents)
        {
            if (a1_.id >= a2_.id || a2_.path.empty())
                continue;

            const auto& a1 = a1_.path.size() <= a2_.path.size() ? a1_ : a2_;
            const auto& a2 = a1_.path.size() <= a2_.path.size() ? a2_ : a1_;
            int t = 1;
            for (; t < (int)a1.path.size(); t++)
            {
                if (a1.path[t].location == a2.path[t].location) // vertex conflict
                {
                    throw ValidationException("Vertex conflict between agents " +
                                              std::to_string(a1.id) + " and " + std::to_string(a2.id) +
                                              " at location " + std::to_string(a1.path[t].location) +
                                              " at timestep " + std::to_string(t));
                }
                else if (a1.path[t].location == a2.path[t - 1].location &&
                         a1.path[t - 1].location == a2.path[t].location) // edge conflict
                {
                    throw ValidationException("Edge conflict between agents " +
                                              std::to_string(a1.id) + " and " + std::to_string(a2.id) +
                                              " at edge (" + std::to_string(a1.path[t - 1].location) + "," +
                                              std::to_string(a1.path[t].location) + ") at timestep " +
                                              std::to_string(t));
                }
            }

            int target = a1.path.back().location;
            for (; t < (int)a2.path.size(); t++)
            {
                if (a2.path[t].location == target)  // target conflict
                {
                    throw ValidationException("Target conflict: agent " + std::to_string(a2.id) +
                                              " (len " + std::to_string(a2.path.size() - 1) + ") traverses target of agent " +
                                              std::to_string(a1.id) + " (len " + std::to_string(a1.path.size() - 1) +
                                              ") at location " + std::to_string(target) + " at timestep " + std::to_string(t));
                }
            }
        }
    }

    if (sum_of_costs != sum)
    {
        throw ValidationException("Sum of costs mismatch: sum_of_costs = " + std::to_string(sum_of_costs) +
                                  ", but computed sum = " + std::to_string(sum));
    }
}

/*
void LNS::validateSolution() const
{
    int sum = 0;
    for (const auto& a1_ : agents)
    {
        if (a1_.path.empty())
        {
            cerr << "No solution for agent " << a1_.id << endl;
            exit(-1);
        }
        else if (a1_.path_planner->start_location != a1_.path.front().location)
        {
            cerr << "The path of agent " << a1_.id << " starts from location " << a1_.path.front().location
                 << ", which is different from its start location " << a1_.path_planner->start_location << endl;
            exit(-1);
        }
        else if (a1_.path_planner->goal_location != a1_.path.back().location)
        {
            cerr << "The path of agent " << a1_.id << " ends at location " << a1_.path.back().location
                 << ", which is different from its goal location " << a1_.path_planner->goal_location << endl;
            exit(-1);
        }
        for (int t = 1; t < (int) a1_.path.size(); t++ )
        {
            if (!instance.validMove(a1_.path[t - 1].location, a1_.path[t].location))
            {
                cerr << "The path of agent " << a1_.id << " jump from "
                     << a1_.path[t - 1].location << " to " << a1_.path[t].location
                     << " between timesteps " << t - 1 << " and " << t << endl;
                exit(-1);
            }
        }
        sum += (int) a1_.path.size() - 1;
        for (const auto  & a2_: agents)
        {
            if (a1_.id >= a2_.id || a2_.path.empty())
                continue;
            const auto & a1 = a1_.path.size() <= a2_.path.size()? a1_ : a2_;
            const auto & a2 = a1_.path.size() <= a2_.path.size()? a2_ : a1_;
            int t = 1;
            for (; t < (int) a1.path.size(); t++)
            {
                if (a1.path[t].location == a2.path[t].location) // vertex conflict
                {
                    cerr << "Find a vertex conflict between agents " << a1.id << " and " << a2.id <<
                         " at location " << a1.path[t].location << " at timestep " << t << endl;
                    exit(-1);
                }
                else if (a1.path[t].location == a2.path[t - 1].location &&
                         a1.path[t - 1].location == a2.path[t].location) // edge conflict
                {
                    cerr << "Find an edge conflict between agents " << a1.id << " and " << a2.id <<
                         " at edge (" << a1.path[t - 1].location << "," << a1.path[t].location <<
                         ") at timestep " << t << endl;
                    exit(-1);
                }
            }
            int target = a1.path.back().location;
            for (; t < (int) a2.path.size(); t++)
            {
                if (a2.path[t].location == target)  // target conflict
                {
                    cerr << "Find a target conflict where agent " << a2.id << " (of length " << a2.path.size() - 1<<
                         ") traverses agent " << a1.id << " (of length " << a1.path.size() - 1<<
                         ")'s target location " << target << " at timestep " << t << endl;
                    exit(-1);
                }
            }
        }
    }
    if (sum_of_costs != sum)
    {
        cerr << "The computed sum of costs " << sum_of_costs <<
             " is different from the sum of the paths in the solution " << sum << endl;
        exit(-1);
    }
}*/



void LNS::writeIterStatsToFile(const string & file_name) const
{
    if (init_lns != nullptr)
    {
        init_lns->writeIterStatsToFile(file_name + "-initLNS.csv");
    }
    if (iteration_stats.size() <= 1)
        return;
    string name = file_name;
    if (use_init_lns or num_of_iterations > 0)
        name += "-LNS.csv";
    else
        name += "-" + init_algo_name + ".csv";
    std::ofstream output;
    output.open(name);
    // header
    output << "num of agents," <<
           "sum of costs," <<
           "runtime," <<
           "cost lowerbound," <<
           "sum of distances," <<
           "MAPF algorithm" << endl;

    for (const auto &data : iteration_stats)
    {
        output << data.num_of_agents << "," <<
               data.sum_of_costs << "," <<
               data.runtime << "," <<
               max(sum_of_costs_lowerbound, sum_of_distances) << "," <<
               sum_of_distances << "," <<
               data.algorithm << endl;
    }
    output.close();
}

void LNS::writeResultToFile(const string & file_name) const
{
    if (init_lns != nullptr)
    {
        init_lns->writeResultToFile(file_name + "-initLNS.csv", sum_of_distances, preprocessing_time);
    }
    string name = file_name;
    if (use_init_lns or num_of_iterations > 0)
        name += "-LNS.csv";
    else
        name += "-" + init_algo_name + ".csv";
    std::ifstream infile(name);
    bool exist = infile.good();
    infile.close();
    if (!exist)
    {
        ofstream addHeads(name);
        addHeads << "runtime,solution cost,initial solution cost,lower bound,sum of distance," <<
                 "iterations," <<
                 "group size," <<
                 "runtime of initial solution,restart times,area under curve," <<
                 "LL expanded nodes,LL generated,LL reopened,LL runs," <<
                 "preprocessing runtime,solver name,instance name" << endl;
        addHeads.close();
    }
    uint64_t num_LL_expanded = 0, num_LL_generated = 0, num_LL_reopened = 0, num_LL_runs = 0;
    for (auto & agent : agents)
    {
        agent.path_planner->reset();
        num_LL_expanded += agent.path_planner->accumulated_num_expanded;
        num_LL_generated += agent.path_planner->accumulated_num_generated;
        num_LL_reopened += agent.path_planner->accumulated_num_reopened;
        num_LL_runs += agent.path_planner->num_runs;
    }
    double auc = 0;
    if (!iteration_stats.empty())
    {
        auto prev = iteration_stats.begin();
        auto curr = prev;
        ++curr;
        while (curr != iteration_stats.end() && curr->runtime < time_limit)
        {
            auc += (prev->sum_of_costs - sum_of_distances) * (curr->runtime - prev->runtime);
            prev = curr;
            ++curr;
        }
        auc += (prev->sum_of_costs - sum_of_distances) * (time_limit - prev->runtime);
    }
    ofstream stats(name, std::ios::app);
    stats << runtime << "," << sum_of_costs << "," << initial_sum_of_costs << "," <<
          max(sum_of_distances, sum_of_costs_lowerbound) << "," << sum_of_distances << "," <<
          iteration_stats.size() << "," << average_group_size << "," <<
          initial_solution_runtime << "," << restart_times << "," << auc << "," <<
          num_LL_expanded << "," << num_LL_generated << "," << num_LL_reopened << "," << num_LL_runs << "," <<
          preprocessing_time << "," << getSolverName() << "," << instance.getInstanceName() << endl;
    stats.close();
}

void LNS::writePathsToFile(const string & file_name) const
{
    std::ofstream output;
    output.open(file_name);
    // header
    // output << agents.size() << endl;

    for (const auto &agent : agents)
    {
        output << "Agent " << agent.id << ":";
        for (const auto &state : agent.path)
            output << "(" << instance.getRowCoordinate(state.location) << "," <<
                            instance.getColCoordinate(state.location) << ")->";
        output << endl;
    }
    output.close();
}
