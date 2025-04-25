#include "LNS.h"
#include "ECBS.h"
#include <queue>
#include <memory>

#include "../include/MAPF.hpp"
#include "SATUtils.h"

LNS::LNS(const Instance& instance, double time_limit, const string & init_algo_name, const string & replan_algo_name,
         const string & destory_name, int neighbor_size, int num_of_iterations, bool use_init_lns,
         const string & init_destory_name, bool use_sipp, int screen, PIBTPPS_option pipp_option) :
         BasicLNS(instance, time_limit, neighbor_size, screen),
         init_algo_name(init_algo_name),  replan_algo_name(replan_algo_name),
         num_of_iterations(num_of_iterations), // nastavuje se v argumentu
         use_init_lns(use_init_lns),init_destory_name(init_destory_name),
         path_table(instance.map_size), pipp_option(pipp_option), metric_rng(std::random_device{}()) {
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
// DESTROY fáze: generateNeighborBySAT() – najde submapu, agenty, T_sync atd.
// --------------------------------------------------------
bool LNS::generateNeighborBySAT() {
    cout << "====================" << endl;
    cout << "SAT destroy operator called." << endl;

    auto [key_agent_id, problematic_timestep] = findBestAgentAndTime();//findMostDelayedAgent();//
    if (key_agent_id < 0) {
        cout << "No delayed agent found." << endl;
        ignored_agents_with_timestep.clear();
        return false;
    }
    cout << "[DEBUG] key_agent_id: " << key_agent_id << endl;
    cout << "[DEBUG] key_agent_id délka globální cesty: " << agents[key_agent_id].path.size() << endl;

    int agent_loc = agents[key_agent_id].path[problematic_timestep].location; // globalID of the cell in 1D matrix
    int submap_size = 25;

    // int agent_id, int submap_size, int agent_location
    auto [submap, agents_in_submap] = getSubmapAndAgents(key_agent_id, submap_size, agent_loc, problematic_timestep);

    std::unordered_set<int> submap_set;
    std::unordered_map<int, pair<int, int>> global_to_local;
    SATUtils::initializeSubmapData(submap, submap_set, global_to_local);

    std::vector<vector<int>> map = SATUtils::generateMapRepresentation(submap, agents_in_submap, problematic_timestep, instance, agents);

    std::vector<int> agents_to_replan = SATUtils::getAgentsToReplan(agents_in_submap, submap_set, problematic_timestep, agents);
    if (agents_to_replan.empty()) {
        cout << "[WARN] No agents to replan in submap." << endl;
        return false;
    }

    int T_sync = problematic_timestep; // synchronizace dle nejproblematičtějšího agenta
    // Debug výpis – lze ponechat nebo odstranit
    vector<pair<int,int>> start_positions, goal_positions;
    for (int agent : agents_to_replan) {
        int start_global = -1, goal_global = -1;
        int goal_time = -1;
        if ((size_t)T_sync >= agents[agent].path.size()) {
            cout << "[ERROR] Agent " << agent << " nemá definovanou pozici v čase T_sync!\n";
            continue;
        }
        int loc_at_Tsync = agents[agent].path[T_sync].location;
        if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
            cout << "[WARNING] Agent " << agent << " není v submapě v čase T_sync!\n";
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
            cout << "[WARNING] Agent " << agent << " nemá platnou cílovou pozici!\n";
            continue;
        }
        auto itS = global_to_local.find(start_global);
        auto itG = global_to_local.find(goal_global);
        if (itS == global_to_local.end() || itG == global_to_local.end()) {
            cout << "[ERROR] Chybí převod globálních souřadnic na lokální!\n";
            continue;
        }
        start_positions.push_back(itS->second);
        goal_positions.push_back(itG->second);
        cout << "Agent " << agent << " | Start (globální): " << start_global
             << " → (lokální): (" << itS->second.first << ", " << itS->second.second << ")"
             << " v čase " << T_sync << " | Cíl (globální): " << goal_global
             << " → (lokální): (" << itG->second.first << ", " << itG->second.second << ")"
             << " v čase " << goal_time << endl;
    }

    neighbor.agents = agents_to_replan;
    neighbor.submap = submap;
    neighbor.submap_set = submap_set;
    neighbor.global_to_local = global_to_local;
    neighbor.map = map;
    neighbor.T_sync = T_sync;
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
/*
int LNS::countConflicts(const Agent& ag)
{
    int c = 0;
    for (int t = 0; t < (int)ag.path.size(); ++t)
        if (path_table.hasConflict(ag.id, ag.path[t].location, t))
            ++c;
    return c;
}*/


// Čím větší skóre, tím „zajímavější“ je agent pro další destrukci/replan.
double LNS::agentScore(const Agent& ag) const {
    const auto& st = ag.stats;
    return  component_weights[0] * st.delay_max +
            component_weights[1] * st.conflict_cnt +
            component_weights[2] * st.stretch_ratio +
            component_weights[3] * (current_iter - st.last_replanned);
}

/*
    1.	Najdi agenta s nejvyšším agentScore(ag).
    2.	U toho agenta vezmi jeho delay_max (maximální zpoždění) a najdi všechny časy t, kde skutečně tohle zpoždění nastalo.
    3.	Ze všech těchto časových kroků náhodně vyber jeden, který jsi ještě neignoroval(a) (sledovat se to musí přes ignored_agents_with_timestep set).
    4.	Pokud už jsou všechny zignorované, vyčisti set a opakuj hledání; pak pokud stále nic, vrátíš {-1,-1} → už nic nedělej.
    Výsledkem je pár (key_agent_id, problematic_timestep).
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

    /* ‑‑ vyber timestep z největších delayů ‑‑ */
    vector<int> ts;
    int dmax = agents[best_id].stats.delay_max;

    /* Najdeme všechny časové kroky t, kde je zpoždění rovno dmax, a které ještě nejsou v setu ignored_agents_with_timestep. */
    for (int t = 0; t < (int)agents[best_id].path.size(); ++t)
    {
        int d = agents[best_id].path_planner->getNumOfDelaysAtTimestep(
                path_table,
                agents[best_id].path,
                agents[best_id].path[t].location,
                t);
        if (d == dmax && !ignored_agents_with_timestep.count({best_id, t}))
            ts.push_back(t);
    }

    /*
      Pokud první průchod nevygeneroval žádný t, znamená to, že jsme všechny jeho nejhorší časy už jednou „ignorovali“.
      Vyčistíme celý ignored_agents_with_timestep (ale jen jednou!) a jednou zopakujeme výběr.
      Pokud ani potom nenajdeme žádný t, vrátíme {-1,-1}, tedy žádná akce.
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
        // 1) maximální zpoždění na jakémkoliv kroku cesty
        st.delay_max     = computeMaxDelay(ag);
        // 2) počet neplatných kroků (edge/vertex konfliktů)
        st.conflict_cnt  = countConflicts(ag);
        // 3) relativní prodloužení cesty vůči heuristickému dolnímu odhadu
        st.stretch_ratio = double(ag.path.size() - 1) /
                           ag.path_planner->my_heuristic[ ag.path_planner->start_location ];
        // 4) poslední iterace, kdy byl daný agent přeplánován
        st.last_replanned = iter;
    }
}

/*
   Chceme, aby metrika, která vedla k úspěchu (tady index 0 = delay_max), rostla, ostatní se slabě decay‑ovaly:
   reaction_factor (např. 0.01) určuje, jak silně se váha „odmění“ za dané δ.
   decay_factor (např. 0.01) určuje, jak rychle se ostatní váhy vrací k větší všestrannosti, když nejsou vybírány.
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

int LNS::selectMetricIndex() const {
    // Discrete distribution that picks i with probability component_weights[i].
    std::discrete_distribution<int> dist(
            component_weights.begin(),
            component_weights.end()
    );
    return dist(metric_rng);
}

// --------------------------------------------------------
// REPAIR fáze: runSAT() – zavolá findLocalPaths + solveWithSAT,
//              a upraví cesty agentů + path_table
// --------------------------------------------------------
bool LNS::runSAT()
{
    cout << "[REPAIR] SAT destroy operator called." << endl;
    cout << "[REPAIR] SAT operator – spouštím subproblém NA OPTIMALIZACI." << endl;

    const auto& agents_to_replan = neighbor.agents;
    const auto& submap           = neighbor.submap;
    const auto& submap_set       = neighbor.submap_set;
    const auto& global_to_local  = neighbor.global_to_local;
    auto& map              = neighbor.map;
    int T_sync                   = neighbor.T_sync;

    // Volání pomocí SATUtils
    auto local_paths = SATUtils::findLocalPaths(agents_to_replan, submap, submap_set, global_to_local, T_sync, agents);

    bool success = SATUtils::solveWithSAT(map, local_paths, agents_to_replan, submap, T_sync, agents);

    if (!success) {
        for (auto a : agents_to_replan)
            path_table.insertPath(agents[a].id, agents[a].path); // return of old paths of agents

        updateAllStats(current_iter); // TODO: zkontrolovat
        cout << "[WARN] SAT solver failed to find a valid solution." << endl;
        return false;
    }
    // --------------------------------------------------------
    //  Dodatečná kontrola platnosti cest vrácených SAT solverem
    //  – pokud agent provádí neplatný skok (např. diagonálu),
    //    vrátíme mu původní cestu z neighbor.old_paths.
    // --------------------------------------------------------
    for (size_t idx = 0; idx < agents_to_replan.size(); ++idx) {
        int ag = agents_to_replan[idx];
        bool invalid_move = false;

        // Projdeme celou (již globální) cestu agenta a ověříme každý krok.
        for (size_t t = 1; t < agents[ag].path.size(); ++t) {
            int from = agents[ag].path[t - 1].location;
            int to   = agents[ag].path[t].location;
            if (!instance.validMove(from, to)) { // zahrnuje i překážky / diag. skoky
                invalid_move = true;
                cout << "[ERROR] Agent " << ag
                     << " má neplatný přechod " << from << " -> " << to
                     << " mezi časy " << t - 1 << " a " << t
                     << ". Vracíme původní cestu." << endl;
                break;
            }
        }

        if (invalid_move) // Obnovíme starou cestu, index v old_paths odpovídá pořadí agents_to_replan.
            agents[ag].path = neighbor.old_paths[idx];
    }

    // úspěch – spočítáme novou sum_of_costs
    neighbor.sum_of_costs = 0;
    for (int ag : agents_to_replan)
        neighbor.sum_of_costs += (int)agents[ag].path.size() - 1;

    if (neighbor.sum_of_costs <= neighbor.old_sum_of_costs) {
        // akceptujeme novou cestu
        for (int a : agents_to_replan) {
            /*
            cout << "[DEBUG] Kontrola STARÉ path_table pro agenta " << a << ":\n";
            for (int t = 0; t < (int) agents[a].path.size(); t++)
            {
                int loc = agents[a].path[t].location;
                if (loc >= 0 && loc < (int) path_table.table.size())
                {
                    // zkontrolovat table[loc].size() > t
                    if ((int) path_table.table[loc].size() > t)
                        cout << "  time=" << t << ", loc=" << loc
                             << ", table=" << path_table.table[loc][t] << endl;
                    else
                        cout << "  time=" << t << ", loc=" << loc << " => out of range\n";
                }
            }*/
            //path_table.deletePath(agents[a].id, agents[a].path);
            path_table.insertPath(agents[a].id, agents[a].path);
            cout << "(LNS.cpp) Nová cesta v agents[a].path agenta " << a << ": ";
            for (auto loc : agents[a].path)
                cout << loc.location << ", ";
            cout << endl;

            /*
            cout << "[DEBUG] Kontrola NOVÉ path_table pro agenta " << a << ":\n";
            for (int t = 0; t < (int) agents[a].path.size(); t++)
            {
                int loc = agents[a].path[t].location;
                if (loc >= 0 && loc < (int) path_table.table.size())
                {
                    if ((int) path_table.table[loc].size() > t)
                        cout << "  time=" << t << ", loc=" << loc
                             << ", table=" << path_table.table[loc][t] << endl;
                    else
                        cout << "  time=" << t << ", loc=" << loc << " => out of range\n";
                }
            }*/
        }

        updateAllStats(current_iter);
        if (neighbor.sum_of_costs < neighbor.old_sum_of_costs) {
            // compute delta based on the chosen metric
            // relativní zlepšení (např. delta = 0.18 = 18% úspora
            double delta = double(neighbor.old_sum_of_costs - neighbor.sum_of_costs)
                           / double(neighbor.old_sum_of_costs);
            cout << "[DEBUG] hodnota delta :" << delta << endl;
            // multiplicative update: vybranému indexu (tady 0) se přidá bonus >1, ostatní se decayí
            int metric_index = selectMetricIndex();
            cout << "[DEBUG] rewarding metric index = " << metric_index << endl;
            updateComponentWeights(metric_index, delta);
            cout << "[DEBUG] component_weights = {"
                 << component_weights[0] << ", "
                 << component_weights[1] << ", "
                 << component_weights[2] << ", "
                 << component_weights[3] << "}" << endl;

        }

        return true;
    } else {
        // revert
        cout << "[INFO] New SAT solution is worse, reverting." << endl;
        for (int i = 0; i < (int)neighbor.agents.size(); i++) {
            int a = neighbor.agents[i];
            //path_table.deletePath(agents[a].id, agents[a].path);
            agents[a].path = neighbor.old_paths[i];
            cout << "(LNS.cpp) Stará cesta v agents[a].path agenta " << a << ": ";
            for (auto loc : agents[a].path)
                cout << loc.location << ", ";
            cout << endl;

            path_table.insertPath(agents[a].id, agents[a].path);
            cout << "[DEBUG] Kontrola STARÉ path_table pro agenta " << a << ":\n";
            for (int t = 0; t < (int) agents[a].path.size(); t++)
            {
                int loc = agents[a].path[t].location;
                if (loc >= 0 && loc < (int) path_table.table.size())
                {
                    // zkontrolovat table[loc].size() > t
                    if ((int) path_table.table[loc].size() > t)
                        cout << "  time=" << t << ", loc=" << loc
                             << ", table=" << path_table.table[loc][t] << endl;
                    else
                        cout << "  time=" << t << ", loc=" << loc << " => out of range\n";
                }
            }
        }
        neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        return false;
    }
}

bool LNS::run()
{
    // Otevřeme soubor pro zápis
    std::ofstream out("log.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(out.rdbuf());

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
                                   replan_algo_name, init_destory_name, neighbor_size, screen);
            succ = init_lns->run(false);
            if (succ)
            {
                path_table.reset();
                for (const auto & agent : agents)
                    path_table.insertPath(agent.id, agent.path);
                init_lns->clear();
                initial_sum_of_costs = init_lns->sum_of_costs;
                sum_of_costs = initial_sum_of_costs;
                updateAllStats(0);
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
    runtime = initial_solution_runtime;
    if (!succ) {
        cout << "Failed to find an initial solution in "
             << runtime << " seconds after  " << restart_times << " restarts" << endl;
        return false;
    }

    // ============================================
    // pomocná lambda funkce pro unify opravu konfliktu
    // ============================================
    auto doInitLNSRepair = [&](const string& debug_reason){
        cout << "[DEBUG] Attempting immediate repair via init_lns " << debug_reason << "." << endl;
        init_lns = new InitLNS(instance, agents, time_limit - runtime,
                               replan_algo_name, init_destory_name, neighbor_size, screen);

        cout << "[DEBUG] In LNS, we pass " << agents.size()
             << " agents to init_lns (skip=true). " << endl;

        // ------------------------------------------------------------------
        // Výpis: path_table pro vybrané agenty (např. neighbor.agents)
        // ------------------------------------------------------------------
        cout << "[DEBUG] Obsah path_table pro vybrané agenty (neighbor.agents):" << endl;
        for (int a : neighbor.agents) {
            cout << "  Agent " << a << " => controlling path length=" << agents[a].path.size() << endl;
            for (int t = 0; t < (int)agents[a].path.size(); t++) {
                int loc = agents[a].path[t].location;
                // zkontrolovat, zda loc je validní index
                if (loc < 0 || loc >= (int)path_table.table.size()) {
                    cout << "    [time=" << t << "]: loc=" << loc << " (out of range)" << endl;
                    continue;
                }
                // zkontrolovat, zda path_table.table[loc].size() > t
                if ((int)path_table.table[loc].size() <= t) {
                    cout << "    [time=" << t << ", loc=" << loc << "]: path_table.table[loc].size()="
                         << path_table.table[loc].size() << " => out of range pro t=" << t << endl;
                    continue;
                }
            }
        }

        bool fixed = init_lns->run(true);

        cout << "[DEBUG] init_lns->sum_of_costs po doběhnutí init_lns->run: " << init_lns->sum_of_costs << endl;

        if (fixed) {
            sum_of_costs = init_lns->sum_of_costs;
            //neighbor.old_sum_of_costs = init_lns->sum_of_costs;

            cout << "[DEBUG] sum_of_costs po přiřazení init_lns->run: " << sum_of_costs << endl;
            //path_table.reset(); // TODO: ověřit
            for (const auto &agent : agents)
                path_table.insertPath(agent.id, agent.path);

            init_lns->clear();
        }
        else cout << "[ERROR] Could not repair solution right after SAT." << endl;
    };
    // ============================================

    bool needConflictRepair = false;

    // optimalizace
    while (runtime < time_limit && iteration_stats.size() <= num_of_iterations) {
        current_iter = static_cast<int>(iteration_stats.size());
        cout.flush();
        runtime = ((fsec)(Time::now() - start_time)).count();

        if (ALNS) chooseDestroyHeuristicbyALNS();

        bool opSuccess = false;
        bool SATchosen = false;

        if (destroy_strategy == SAT) { // SAT is destroy operator merged with replan operator
            int r = rand() % 100;
            if (r < 100) { // číslo zde bude hyperparametr
                SATchosen = true;
                //cout << "[DEBUG] hodnota r je " << r << endl;
                cout << "[DEBUG] Using SAT operator (destroy+repair SAT)." << endl;
                //const int MAX_SAT_ATTEMPTS = 10;
                while (!opSuccess) {
                    if (!generateNeighborBySAT()) continue;
                    neighbor.old_paths.resize(neighbor.agents.size());
                    neighbor.old_sum_of_costs = 0;
                    for (int i = 0; i < (int)neighbor.agents.size(); i++) {
                        int a = neighbor.agents[i];
                        neighbor.old_paths[i] = agents[a].path;
                        path_table.deletePath(a, agents[a].path);
                        neighbor.old_sum_of_costs += (int)agents[a].path.size() - 1;
                    }

                    // Removed stats and component weights updates from run() to be handled in runSAT().

                    opSuccess = runSAT();
                }
            }
            else cout << "[DEBUG] Random chance did not select SAT operator (r=" << r << "), using default strategy." << endl;
            cout << "[DEBUG] hodnota opSuccess: " << opSuccess << endl;
        }

        if (!opSuccess)
        {
            // fallback neighbor generation
            int DEFAULT_DESTROY_STRATEGY = INTERSECTION;
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
                    cerr << "Wrong neighbor generation strategy" << endl;
                    exit(-1);
            }

            if (!opSuccess)
                continue;

            neighbor.old_paths.resize(neighbor.agents.size());
            neighbor.old_sum_of_costs = 0;
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                int a = neighbor.agents[i];
                neighbor.old_paths[i] = agents[a].path;
                path_table.deletePath(a, agents[a].path);
                neighbor.old_sum_of_costs += (int)agents[a].path.size() - 1;
            }

            std::string DEFAULT_REPLAN_ALGO = "PP"; // fixně
            if      (DEFAULT_REPLAN_ALGO == "PP")   succ = runPP();
            else if (DEFAULT_REPLAN_ALGO == "CBS")  succ = runCBS();
            else if (DEFAULT_REPLAN_ALGO == "EECBS")succ = runEECBS();
            else { cerr << "Wrong replanning strategy" << endl; exit(-1); }
        }
        else succ = opSuccess;// opSuccess = true => runSAT proběhl

        if (!succ)
            continue;

        // ALNS vyhodnocení
        if (ALNS)
        {
            if (neighbor.old_sum_of_costs > neighbor.sum_of_costs)
                destroy_weights[selected_neighbor] =
                        reaction_factor * (neighbor.old_sum_of_costs - neighbor.sum_of_costs) / neighbor.agents.size()
                        + (1 - reaction_factor) * destroy_weights[selected_neighbor];
            else
                destroy_weights[selected_neighbor] =
                        (1 - decay_factor) * destroy_weights[selected_neighbor];
        }

        runtime = ((fsec)(Time::now() - start_time)).count();

        cout << "[DEBUG] neighbor.sum_of_costs před opětovným přepočtem: " << neighbor.sum_of_costs << endl;
        cout << "[DEBUG] neighbor.old_sum_of_costs před opětovným přepočtem: " << neighbor.old_sum_of_costs << endl;
        cout << "[DEBUG] sum_of_costs před opětovným přepočtem: " << sum_of_costs << endl;

        // ------------------------------------------------
        // 2) Po SAT => validace a případný conflict repair
        // ------------------------------------------------
        if (destroy_strategy == SAT && opSuccess && SATchosen)
        {
            /* ---  synchronize global sum_of_costs *před* validací  ---
               runSAT právě změnil cesty vybraných agentů a naplnil
               neighbor.{old_,}sum_of_costs.  validateSolution() kontroluje
               konzistenci proměnné sum_of_costs, takže ji musíme
               přepočítat dřív, jinak o 1/2 kroky zaostává a hlásí
               “sum of costs mismatch”.                                          */
            sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;
            /* aby se o pár řádků níže delta nepřičetla znovu */
            neighbor.old_sum_of_costs = neighbor.sum_of_costs;

            cout << "[DEBUG] Validate solution immediately after SAT success." << endl;
            try {
                validateSolution();
                cout << "[DEBUG] No problems after SAT replan." << endl;
            } catch (const ValidationException& e) {
                cout << "[WARNING] Problem after SAT: " << e.what() << endl;
                // unify
                doInitLNSRepair("because problem occured after SAT (should be applied only for conflicts...)");
                continue;
            }
        } else sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;

        cout << "[DEBUG] recomputing sum_of_cost by dividing neighbor.sum_of_costs and neighbor.old_sum_of_costs" << endl;
        cout << "[DEBUG] sum_of_costs po opětovném přepočtu: " << sum_of_costs << endl;

        cout << "Iteration " << iteration_stats.size() << endl;

        if (screen >= 1)
        {
            cout << "Iteration " << iteration_stats.size()
                 << ", group size = " << neighbor.agents.size()
                 << ", solution cost = " << sum_of_costs
                 << ", remaining time = " << time_limit - runtime
                 << endl;
        }
        iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, replan_algo_name);
    }

    average_group_size = -iteration_stats.front().num_of_agents;
    for (const auto& data : iteration_stats)
        average_group_size += data.num_of_agents;
    if (average_group_size > 0)
        average_group_size /= (double)(iteration_stats.size() - 1);

    cout << getSolverName()
         << ": runtime = " << runtime
         << ", iterations = " << iteration_stats.size()
         << ", solution cost = " << sum_of_costs
         << ", initial solution cost = " << initial_sum_of_costs
         << ", failed iterations = " << num_of_failures
         << endl;

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
    else if (init_algo_name == "SAT") // TODO: popravdě nevím, jestli máme upravovat - možná není nutné
        runPP();
        //succ = runSAT();
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
        updateAllStats(current_iter); // PŘIDÁNO
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
        case 3 : destroy_strategy = SAT; break;
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

    std::pair<int, int> result = findMostDelayedAgent();
    int a = result.first;
    int time_step = result.second;

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

/* PRO FUNGOVÁNÍ RandomWalk
 * int LNS::findMostDelayedAgent()
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
}*/

pair<int,int> LNS::findMostDelayedAgent()
{
    const int N = static_cast<int>(agents.size());
    if (N == 0) return {-1,-1};

    /* Projdeme agenty cyklicky: (last+1) … (last+N) */
    for (int step = 1; step <= N; ++step)
    {
        int idx = (last_selected_agent + step) % N;
        const auto& ag = agents[idx];

        // najdi NEJvíc zpožděný timestep tohoto agenta,
        // který zatím není v ignored_agents_with_timestep
        auto [delays, ts] =
                ag.getMostProblematicDelay(path_table, ignored_agents_with_timestep);

        if (ts < 0)           // už nemá nic nového
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


/*
pair<int, int> LNS::findMostDelayedAgent() {
    int max_delays = -1;
    int agent_with_max_delays = -1;
    int most_problematic_timestep = -1;

    for (const auto& agent : agents) {
        // přeskoč agenty, kteří už selhali
        // TODO: tohle není nejmoudřejší, protože je nechceme přeskakovat natrvalo, ale prozatím to nevadí
        // TODO: nikdy nevymazáváme ignored_agents.clear()
        auto [agent_max_delays, problematic_timestep] = // TODO: zvolit ještě jiný způsob výběru
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
}*/

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
