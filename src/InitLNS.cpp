#include "InitLNS.h"
#include <queue>
#include <algorithm>
#include "GCBS.h"
#include "PBS.h"

#include "../include/MAPF.hpp"
#include "SATUtils.h"

InitLNS::InitLNS(const Instance& instance, vector<Agent>& agents, double time_limit,
         const string & replan_algo_name, const string & init_destory_name, int neighbor_size, int screen) :
         BasicLNS(instance, time_limit, neighbor_size, screen), agents(agents), replan_algo_name(replan_algo_name),
         path_table(instance.map_size, agents.size()), collision_graph(agents.size()), goal_table(instance.map_size, -1) {
     replan_time_limit = time_limit;
     if (init_destory_name == "Adaptive")
     {
         ALNS = true;
         destroy_weights.assign(INIT_COUNT, 1);
         decay_factor = 0.05;
         reaction_factor = 0.05;
     }
     else if (init_destory_name == "Target")
         init_destroy_strategy = TARGET_BASED;
     else if (init_destory_name == "Collision")
         init_destroy_strategy = COLLISION_BASED;
     else if (init_destory_name == "Random")
         init_destroy_strategy = RANDOM_BASED;
     else if (init_destory_name == "SAT")
         init_destroy_strategy = SAT_BASED;
     else
     {
         cerr << "Init Destroy heuristic " << init_destory_name << " does not exists. " << endl;
         exit(-1);
     }

     for (auto& i:agents) {
         goal_table[i.path_planner->goal_location] = i.id;
     }
}

/*
* findConflictAgent
*  - vrací agenta a časový krok, kdy je poprvé zjištěn konflikt
*  - Pokud collision_graph[a].size() == 0 pro všechny a, vrátí {-1,-1} => žádné konflikty.
*  - Pro zjednodušení vracíme conflict_time=0;
*    v reálné implementaci ho lze najít analýzou path_tableWC (např. hasCollisions / getLastCollisionTimestep).
*/
pair<int, int> InitLNS::findConflictAgent() {
    cout << "[DEBUG] Agenti v sat_failed_agents: ";
    for (auto a : failed_sat_agents)
        cout << "agent " << a << ", ";
    cout << endl;
    for (const auto& agent : agents) {
        if (failed_sat_agents.find(agent.id) != failed_sat_agents.end())
            continue;
        if (collision_graph[agent.id].empty())
            continue;
        const auto& path = agent.path;
        if (path.empty())
            continue;

        for (int t = 1; t < (int)path.size(); t++) { // t=1 kvůli edge konfliktům
            int from = path[t - 1].location;
            int to = path[t].location;

            if (path_table.hasCollisions(from, to, t))
                return {agent.id, t-1}; // agent a čas konfliktu - 1 (kvůli vertex konfliktům)
        }
    }

    return {-1, -1}; // žádný konflikt nenalezen
}

// getting the submap around one or more agents and identifying agents in these submaps
pair<vector<vector<int>>, vector<int>> InitLNS::getSubmapAndAgents(int agent_id, int submap_size, int agent_location, int timestep) {
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
    int agent_x = agent_location / map_width;
    int agent_y = agent_location % map_width;

    int half_side = submap_side / 2;

    for (int dx = -half_side; dx <= half_side; ++dx) {
        for (int dy = -half_side; dy <= half_side; ++dy) {
            int x = agent_x + dx;
            int y = agent_y + dy;

            // ensure we are within map boundaries
            if (x >= 0 && x < map_height && y >= 0 && y < map_width) {
                int global_pos = x * map_width + y; // unique index of cell in global map

                // map (dx, dy) to submap indices
                int submap_x = dx + half_side;
                int submap_y = dy + half_side;

                if (submap_x >= 0 && submap_x < submap_side && submap_y >= 0 && submap_y < submap_side) {
                    submap[submap_x][submap_y] = global_pos;
                    path_table.get_agents_at_timestep(conflicting_agents, global_pos, timestep);
                    //path_table.get_agents(conflicting_agents, global_pos);
                }
            }
        }
    }

    agents_in_submap.assign(conflicting_agents.begin(), conflicting_agents.end());
    bool found = false;
    for (auto a : agents_in_submap) {
        if (a == agent_id) {
            found = true;
            break;
        }
    }
    if (found)
        cout << "[DEBUG] klíčový agent je mezi agenty v submapě" << endl;
    else cout << "[WARNING] klíčový agent NENÍ mezi agenty v submapě!" << endl;
    return {submap, agents_in_submap};
}

// --------------------------------------------------------
// DESTROY fáze: generateNeighborBySAT() – najde submapu, agenty, T_sync atd.
// --------------------------------------------------------
bool InitLNS::generateNeighborBySAT() {
    cout << "====================" << endl;
    cout << "SAT destroy operator called." << endl;

    auto [key_agent_id, problematic_timestep] = findConflictAgent();

    //problematic_timestep += 10; // DOČASNĚ PRO TEST

    if (key_agent_id < 0) {
        cout << "[DEBUG] Žádný agent s konflikty nebyl nalezen." << endl;
        return false;
    }
    else {
        cout << "[DEBUG] Vybraný agent " << key_agent_id << " má "
             << collision_graph[key_agent_id].size() << " konfliktů." << endl;
    }

    cout << "[DEBUG] --- Speciální debug pro key_agent_id ---" << endl;
    cout << "[DEBUG] Collision graph pro agenta key_agent_id: ";
    for (int conflict : collision_graph[key_agent_id])
        cout << conflict << " ";
    cout << endl;
    cout << "[DEBUG] Problematický timestep: " << problematic_timestep << endl;

    int agent_loc = agents[key_agent_id].path[problematic_timestep].location; // globalID buňky
    int submap_size = 9;

    // Získání submapy a seznamu agentů
    auto [submap, agents_in_submap] = getSubmapAndAgents(key_agent_id, submap_size, agent_loc, problematic_timestep);

    cout << "[DEBUG] Počet agentů v submapě: " << agents_in_submap.size() << endl;
    if (!agents_in_submap.empty()) {
        cout << "[DEBUG] Seznam agentů v submapě: ";
        for (int ag : agents_in_submap)
            cout << ag << " ";
        cout << endl;
    }
    else {
        cout << "[DEBUG] Upozornění: Submapa je prázdná, žádní agenti nebyli nalezeni v okolí." << endl;
    }

    // Připravíme data pro SATUtils
    std::unordered_set<int> submap_set;
    std::unordered_map<int, pair<int, int>> global_to_local;
    SATUtils::initializeSubmapData(submap, submap_set, global_to_local);

    // Debug – výpis obsahu submap_set
    cout << "[DEBUG] Obsah submap_set (globální indexy): ";
    for (int pos : submap_set)
        cout << pos << " ";
    cout << endl;

    // Generování 2D reprezentace mapy s překážkami a umístěním agentů
    vector<vector<int>> map = SATUtils::generateMapRepresentation(submap, agents_in_submap, problematic_timestep, instance, agents);

    cout << "[DEBUG] Kontrola: Agent " << key_agent_id << " má pozici " << agent_loc
         << " a T_sync = " << problematic_timestep << endl;

    // Volání SATUtils::getAgentsToReplan
    vector<int> agents_to_replan = SATUtils::getAgentsToReplan(agents_in_submap, submap_set, problematic_timestep, agents);
    if (agents_to_replan.empty()) {
        cout << "[WARN] No agents to replan in submap." << endl;
        cout << "[DEBUG] Pro agenta key_agent_id: agents_in_submap je prázdný. Kontrola submap_set a problematic timestep:" << endl;
        cout << "         T_sync: " << problematic_timestep << ", agent_loc: " << agent_loc << endl;
        return false;
    }

    int T_sync = problematic_timestep; // synchronizace dle nejkonfliktnějšího agenta

    // Debug – výpis start/goal pozic pro agenty k přeplánování
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
            cout << "[WARNING] Agent " << agent << " nemá platnou cílovou pozici (už T_sync je mimo?)!\n";
            continue;
        }
        auto itS = global_to_local.find(start_global);
        auto itG = global_to_local.find(goal_global);
        if (itS == global_to_local.end() || itG == global_to_local.end()) {
            cout << "[ERROR] Převod globálních souřadnic selhal pro agenta " << agent << "!\n";
            continue;
        }
        start_positions.push_back(itS->second);
        goal_positions.push_back(itG->second);
        cout << "Agent " << agent
             << " | Start (globální): " << start_global
             << " → (lokálně): (" << itS->second.first << ", " << itS->second.second << ")"
             << " v čase " << T_sync
             << " | První opuštění submapy v t=" << (goal_time + 1)
             << " => Cíl (globální): " << goal_global
             << " → (lokálně): (" << itG->second.first << ", " << itG->second.second << ")"
             << " v čase " << goal_time << endl;
    }

    // Uložení dat do neighbor pro pozdější repair fázi
    neighbor.agents = agents_to_replan;
    neighbor.submap = submap;
    neighbor.submap_set = submap_set;
    neighbor.global_to_local = global_to_local;
    neighbor.map = map;
    neighbor.T_sync = T_sync;
    neighbor.key_agent_id = key_agent_id;

    return true;
}

// TODO: špatně se vybírá konfliktní agent nebo se špatně generuje submapa

// --------------------------------------------------------
// REPAIR fáze: runSAT() – zavolá findLocalPaths + solveWithSAT,
//              a upraví cesty agentů + path_table
// --------------------------------------------------------
bool InitLNS::runSAT()
{
    cout << "====================" << endl;
    cout << "[REPAIR] SAT operator – spouštím subproblém NA ŘEŠENÍ KONFLIKTŮ." << endl;

    const auto& agents_to_replan = neighbor.agents;
    const auto& submap = neighbor.submap;
    const auto& submap_set = neighbor.submap_set;
    const auto& global_to_local = neighbor.global_to_local;
    auto& map = neighbor.map;
    int T_sync = neighbor.T_sync;
    int key_agent_id = neighbor.key_agent_id;

    // Voláme SATUtils::findLocalPaths – vytvoří lokální cesty z globálních pomocí předaných dat
    auto local_paths = SATUtils::findLocalPaths(agents_to_replan, submap, submap_set, global_to_local, T_sync, agents);

    // Voláme SATUtils::solveWithSAT – přeplánujeme lokální cesty
    bool success = SATUtils::solveWithSAT(map, local_paths, agents_to_replan, submap, T_sync, agents);

    if (!success) {
        cout << "[WARN] SAT solver failed to find a valid solution." << endl;
        // Vrátíme původní cesty agentům, aby jejich přeplánování nezmizelo.
        for (int i = 0; i < (int)neighbor.agents.size(); i++) {
            int a = neighbor.agents[i];
            cout << "[DEBUG] Reverting path for agent " << a << " (agent id: " << agents[a].id << ")" << endl;
            // Odstraníme aktuální (neúspěšnou) cestu z path_table
            //path_table.deletePath(agents[a].id);
            // Obnovíme původní cestu uloženou v neighbor.old_paths
            agents[a].path = neighbor.old_paths[i];
            // Znovu vložíme původní cestu do path_table
            path_table.insertPath(agents[a].id, agents[a].path);
        }
        failed_sat_agents.insert(key_agent_id);
        // Můžete také vymazat key_agent_id z failed_sat_agents nebo jej nepřidávat
        return false;
    }

    set<pair<int,int>> new_colliding_pairs;
    for (int ag : agents_to_replan)
        updateCollidingPairs(new_colliding_pairs, agents[ag].id, agents[ag].path);
    int new_conflicts = new_colliding_pairs.size();
    cout << "[DEBUG] Nově přeplánované řešení má " << new_conflicts << " konfliktů." << endl;
    cout << "[DEBUG] Původní počet konfliktů: " << neighbor.old_colliding_pairs.size() << endl;

    // porovnáme s původním počtem kolizních párů, uloženým v neighbor.old_colliding_pairs
    if (new_conflicts <= neighbor.old_colliding_pairs.size()) {
        // Akceptujeme nové řešení – aktualizujeme path_table
        for (int a : agents_to_replan) {
            cout << "[DEBUG] we decided to accept the solution -> replanning" << endl;
            cout << "[DEBUG] before inserting new path into path_table, path_table has agent " << a
                 << " with path length = " << path_table.getPath(agents[a].id)->size() << endl;
            path_table.insertPath(agents[a].id, agents[a].path);
            cout << "[DEBUG] after inserting new path into path_table, agents " << agents[a].id
                 << " path has path length = " << path_table.getPath(agents[a].id)->size() << endl;
            cout << "[DEBUG] after inserting new path into path_table, agents " << agents[a].id
                 << " path has path length in agents[a].path = " << agents[a].path.size() << endl;
        }
        failed_sat_agents.clear();
        return true;
    }
    else {
        cout << "[INFO] New SAT solution has more conflicts (" << new_conflicts
             << ") than before (" << neighbor.old_colliding_pairs.size() << "), reverting." << endl;
        // Vrátíme staré cesty
        for (int i = 0; i < (int)neighbor.agents.size(); i++) {
            int a = neighbor.agents[i];
            //path_table.deletePath(agents[a].id);
            agents[a].path = neighbor.old_paths[i];
            path_table.insertPath(agents[a].id, agents[a].path);
            auto storedPath = path_table.getPath(agents[a].id);
            if (storedPath == nullptr)
                cout << "[ERROR] after revert, path_table says agent " << a << " is nullptr" << endl;
            else {
                cout << "[DEBUG] after revert, path_table has agent " << a
                     << " with path length = " << storedPath->size() << endl;
                cout << "[DEBUG] after revert, agents " << agents[a].id
                        << " path has path length = " << agents[a].path.size() << endl;
            }

        }
        failed_sat_agents.insert(key_agent_id);
        return false;
    }
}

bool InitLNS::run()
{
    start_time = Time::now();
    bool succ = getInitialSolution();
    runtime = ((fsec)(Time::now() - start_time)).count();
    iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, "PP", 0, num_of_colliding_pairs);
    if (screen >= 3)
        printPath();
    if (screen >= 1)
        cout << "Iteration " << iteration_stats.size() << ", "
             << "group size = " << neighbor.agents.size() << ", "
             << "colliding pairs = " << num_of_colliding_pairs << ", "
             << "solution cost = " << sum_of_costs << ", "
             << "remaining time = " << time_limit - runtime << endl;
    if (runtime >= time_limit && !succ)
    {
        printResult();
        return false;
    }

    vector<Path*> paths(agents.size());
    for (size_t i = 0; i < agents.size(); i++)
        paths[i] = &agents[i].path;

    while (runtime < time_limit && num_of_colliding_pairs > 0)
    {
        assert(instance.validateSolution(paths, sum_of_costs, num_of_colliding_pairs));
        if (ALNS)
            chooseDestroyHeuristicbyALNS();

        // Náhodně vybereme operátor: s 20% pravděpodobností použijeme SAT operátor,
        // jinak použijeme ostatní operátory dle strategie (TARGET, COLLISION, RANDOM).
        int r = rand() % 100;
        bool opSuccess = false;
        if (r < 20) {
            cout << "[DEBUG] Using SAT operator (destroy+repair SAT) with probability " << r << " %." << endl;
            const int MAX_SAT_ATTEMPTS = 10;
            bool sat_success = false;
            for (int attempt = 0; attempt < MAX_SAT_ATTEMPTS && !sat_success; attempt++) {
                if (!generateNeighborBySAT())
                    continue; // nepodařilo se nalézt validní neighborhood – zkuste znovu

                // najdi kolizní páry
                neighbor.old_colliding_pairs.clear();
                for (int a : neighbor.agents)
                    for (auto j: collision_graph[a])
                        neighbor.old_colliding_pairs.emplace(min(a, j), max(a, j));

                // uložení starých cest
                neighbor.old_paths.resize(neighbor.agents.size());
                neighbor.old_sum_of_costs = 0;
                for (int i = 0; i < (int)neighbor.agents.size(); i++) {
                    int a = neighbor.agents[i];
                    neighbor.old_paths[i] = agents[a].path;
                    // ============ DŮLEŽITÉ ============
                    // tady smažeme agentům cestu z path_table, buď ji po přeplánování nahradíme novou, nebo vrátíme starou
                    // tohle nahrazování nových cest v path_table děláme na konci v runSAT
                    path_table.deletePath(agents[a].id);
                    neighbor.old_sum_of_costs += (int) agents[a].path.size() - 1;
                }

                if (runSAT())
                    sat_success = true;
            }
            opSuccess = sat_success;
        } else {
            cout << "[DEBUG] zavolalo se něco jiného, konkrétně ";
            // Pokud byla zvolena SAT strategie už uživatelem, nahrazujeme ji náhodným výběrem z TARGET, COLLISION, RANDOM
            int strategy;
            if (init_destroy_strategy == SAT_BASED)
                strategy = rand() % 3; // 0: TARGET_BASED, 1: COLLISION_BASED, 2: RANDOM_BASED
            else
                strategy = init_destroy_strategy;

            cout << strategy << endl;
            switch (strategy)
            {
                case TARGET_BASED:
                    opSuccess = generateNeighborByTarget();
                    break;
                case COLLISION_BASED:
                    opSuccess = generateNeighborByCollisionGraph();
                    break;
                case RANDOM_BASED:
                    opSuccess = generateNeighborRandomly();
                    break;
                default:
                    cerr << "Wrong neighbor generation strategy" << endl;
                    exit(-1);
            }
            if (!opSuccess)
                continue;

            // najdi kolizní páry
            neighbor.old_colliding_pairs.clear();
            for (int a : neighbor.agents)
                for (auto j: collision_graph[a])
                    neighbor.old_colliding_pairs.emplace(min(a, j), max(a, j));

            // uložení starých cest
            neighbor.old_paths.resize(neighbor.agents.size());
            neighbor.old_sum_of_costs = 0;
            for (int i = 0; i < (int)neighbor.agents.size(); i++) {
                int a = neighbor.agents[i];
                neighbor.old_paths[i] = agents[a].path;
                // ============ DŮLEŽITÉ ============
                // tady smažeme agentům cestu z path_table, buď ji po přeplánování nahradíme novou, nebo vrátíme starou
                // tohle nahrazování nových cest v path_table děláme na konci v runSAT
                path_table.deletePath(agents[a].id);
                neighbor.old_sum_of_costs += (int) agents[a].path.size() - 1;
            }
        }

        if (!opSuccess || neighbor.agents.empty())
            continue;

        // Replan podle zvolené strategie
        if (replan_algo_name == "PP" || neighbor.agents.size() == 1)
            succ = runPP();
        else if (replan_algo_name == "GCBS")
            succ = runGCBS();
        else if (replan_algo_name == "SAT") {
            // SAT replanning jsme již provedli výše
        }
        else if (replan_algo_name == "PBS")
            succ = runPBS();
        else {
            cerr << "Wrong replanning strategy" << endl;
            exit(-1);
        }

        if (ALNS)
        {
            if (neighbor.colliding_pairs.size() < neighbor.old_colliding_pairs.size())
                destroy_weights[selected_neighbor] =
                        reaction_factor * (double)(neighbor.old_colliding_pairs.size() - neighbor.colliding_pairs.size())
                        + (1 - reaction_factor) * destroy_weights[selected_neighbor];
            else
                destroy_weights[selected_neighbor] =
                        (1 - decay_factor) * destroy_weights[selected_neighbor];
        }

        if (screen >= 2)
            cout << "New colliding pairs = " << neighbor.colliding_pairs.size() << endl;
        if (succ)
        {
            num_of_colliding_pairs += (int)neighbor.colliding_pairs.size() - (int)neighbor.old_colliding_pairs.size();
            for (const auto& agent_pair : neighbor.old_colliding_pairs)
            {
                collision_graph[agent_pair.first].erase(agent_pair.second);
                collision_graph[agent_pair.second].erase(agent_pair.first);
            }
            for (const auto& agent_pair : neighbor.colliding_pairs)
            {
                collision_graph[agent_pair.first].emplace(agent_pair.second);
                collision_graph[agent_pair.second].emplace(agent_pair.first);
            }
            if (screen >= 2)
                printCollisionGraph();
        }

        runtime = ((fsec)(Time::now() - start_time)).count();
        sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;
        if (screen >= 1)
            cout << "Iteration " << iteration_stats.size() << ", "
                 << "group size = " << neighbor.agents.size() << ", "
                 << "colliding pairs = " << num_of_colliding_pairs << ", "
                 << "solution cost = " << sum_of_costs << ", "
                 << "remaining time = " << time_limit - runtime << endl;
        iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, replan_algo_name,
                                     0, num_of_colliding_pairs);
    }

    printResult();
    return (num_of_colliding_pairs == 0);
}

bool InitLNS::runGCBS()
{
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(agents[i].path_planner);
    }

    // build path tables
    vector<PathTable> path_tables(neighbor.agents.size(), PathTable(instance.map_size));
    for (int i = 0; i < (int)neighbor.agents.size(); i++)
    {
        int agent_id = neighbor.agents[i];
        for (int j = 0; j < instance.getDefaultNumberOfAgents(); j++)
        {
            if (j != agent_id and collision_graph[agent_id].count(j) == 0)
                path_tables[i].insertPath(j, agents[j].path);
        }
    }

    GCBS gcbs(search_engines, screen - 1, &path_tables);
    gcbs.setDisjointSplitting(false);
    gcbs.setBypass(true);
    gcbs.setTargetReasoning(true);

    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime;
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    gcbs.solve(T);
    if (gcbs.best_node->colliding_pairs < (int) neighbor.old_colliding_pairs.size()) // accept new paths
    {
        auto id = neighbor.agents.begin();
        neighbor.colliding_pairs.clear();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *gcbs.paths[i];
            updateCollidingPairs(neighbor.colliding_pairs, agents[*id].id, agents[*id].path);
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = gcbs.best_node->sum_of_costs;
        return true;
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
        num_of_failures++;
        return false;
    }
}
bool InitLNS::runPBS()
{
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    vector<const Path*> initial_paths;
    initial_paths.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(agents[i].path_planner);
        initial_paths.push_back(&agents[i].path);
    }

    PBS pbs(search_engines, path_table, screen - 1);
    // pbs.setInitialPath(initial_paths);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime;
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = pbs.solve(T, (int)neighbor.agents.size(), neighbor.old_colliding_pairs.size());
    if (succ and pbs.best_node->getCollidingPairs() < (int) neighbor.old_colliding_pairs.size()) // accept new paths
    {
        auto id = neighbor.agents.begin();
        neighbor.colliding_pairs.clear();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *pbs.paths[i];
            updateCollidingPairs(neighbor.colliding_pairs, agents[*id].id, agents[*id].path);
            path_table.insertPath(agents[*id].id);
            ++id;
        }
        assert(neighbor.colliding_pairs.size() == pbs.best_node->getCollidingPairs());
        neighbor.sum_of_costs = pbs.best_node->sum_of_costs;
        return true;
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id);
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
        num_of_failures++;
        return false;
    }
}
bool InitLNS::runPP()
{
    cout << "[DEBUG] velikost neighbor.agents: " << neighbor.agents.size() << endl;
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());
    if (screen >= 2) {
        cout<<"Neighbors_set: ";
        for (auto id : shuffled_agents)
            cout << id << ", ";
        cout << endl;
    }
    int remaining_agents = (int)shuffled_agents.size();
    auto p = shuffled_agents.begin();
    neighbor.sum_of_costs = 0;
    neighbor.colliding_pairs.clear();
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = min(time_limit - runtime, replan_time_limit);
    auto time = Time::now();
    ConstraintTable constraint_table(instance.num_of_cols, instance.map_size, nullptr, &path_table);
    while (p != shuffled_agents.end() && ((fsec)(Time::now() - time)).count() < T)
    {
        if (p == shuffled_agents.end()) {
            cerr << "[ERROR] Iterator p je na konci shuffled_agents. shuffled_agents.size() = "
                 << shuffled_agents.size() << endl;
            exit(-1);
        }
        int id = *p;
        if (agents[id].path_planner == nullptr) {
            cerr << "[ERROR] agents[" << id << "].path_planner je nullptr!" << endl;
            exit(-1); // nebo vhodně ošetřit chybu
        }
        agents[id].path = agents[id].path_planner->findPath(constraint_table);
        assert(!agents[id].path.empty() && agents[id].path.back().location == agents[id].path_planner->goal_location);
        if (agents[id].path_planner->num_collisions > 0)
            updateCollidingPairs(neighbor.colliding_pairs, agents[id].id, agents[id].path);

        assert(agents[id].path_planner->num_collisions > 0 or
            !updateCollidingPairs(neighbor.colliding_pairs, agents[id].id, agents[id].path));
        neighbor.sum_of_costs += (int)agents[id].path.size() - 1;
        remaining_agents--;

        if (screen >= 3)
        {
            runtime = ((fsec)(Time::now() - start_time)).count();
            cout << "After agent " << id << ": Remaining agents = " << remaining_agents <<
                 ", colliding pairs = " << neighbor.colliding_pairs.size() <<
                 ", LL nodes = " << agents[id].path_planner->getNumExpanded() <<
                 ", remaining time = " << time_limit - runtime << " seconds. " << endl;
        }
        if (neighbor.colliding_pairs.size() > neighbor.old_colliding_pairs.size())
            break;
        path_table.insertPath(agents[id].id, agents[id].path);
        ++p;
    }
    if (p == shuffled_agents.end() && neighbor.colliding_pairs.size() <= neighbor.old_colliding_pairs.size()) // accept new paths
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
            path_table.deletePath(agents[a].id);
            ++p2;
        }
        if (!neighbor.old_paths.empty())
        {
            p2 = neighbor.agents.begin();
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                int a = *p2;
                agents[a].path = neighbor.old_paths[i];
                path_table.insertPath(agents[a].id);
                ++p2;
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
        return false;
    }
}

bool InitLNS::getInitialSolution()
{
    neighbor.agents.clear();
    neighbor.agents.reserve(agents.size());
    sum_of_costs = 0;
    for (int i = 0; i < (int)agents.size(); i++)
    {
        if (agents[i].path.empty())
            neighbor.agents.push_back(i);
        else
        {
            sum_of_costs += (int)agents[i].path.size() - 1;
            path_table.insertPath(agents[i].id, agents[i].path);
        }
    }
    int remaining_agents = (int)neighbor.agents.size();
    std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
    ConstraintTable constraint_table(instance.num_of_cols, instance.map_size, nullptr, &path_table);
    set<pair<int, int>> colliding_pairs;
    for (auto id : neighbor.agents)
    {
        agents[id].path = agents[id].path_planner->findPath(constraint_table);
        assert(!agents[id].path.empty() && agents[id].path.back().location == agents[id].path_planner->goal_location);
        if (agents[id].path_planner->num_collisions > 0)
            updateCollidingPairs(colliding_pairs, agents[id].id, agents[id].path);
        sum_of_costs += (int)agents[id].path.size() - 1;
        remaining_agents--;
        path_table.insertPath(agents[id].id, agents[id].path);
        runtime = ((fsec)(Time::now() - start_time)).count();
        if (screen >= 3)
        {
            cout << "After agent " << id << ": Remaining agents = " << remaining_agents <<
                 ", colliding pairs = " << colliding_pairs.size() <<
                 ", LL nodes = " << agents[id].path_planner->getNumExpanded() <<
                 ", remaining time = " << time_limit - runtime << " seconds. " << endl;
        }
        if (runtime > time_limit)
            break;
    }

    num_of_colliding_pairs = colliding_pairs.size();
    for(const auto& agent_pair : colliding_pairs)
    {
        collision_graph[agent_pair.first].emplace(agent_pair.second);
        collision_graph[agent_pair.second].emplace(agent_pair.first);
    }
    if (screen >= 2)
        printCollisionGraph();
    return remaining_agents == 0;
}

// return true if the new p[ath has collisions;
bool InitLNS::updateCollidingPairs(set<pair<int, int>>& colliding_pairs, int agent_id, const Path& path) const
{
    bool succ = false;
    if (path.size() < 2)
        return succ;
    for (int t = 1; t < (int)path.size(); t++)
    {
        int from = path[t - 1].location;
        int to = path[t].location;
        if ((int)path_table.table[to].size() > t) // vertex conflicts
        {
            for (auto id : path_table.table[to][t])
            {
                succ = true;
                colliding_pairs.emplace(min(agent_id, id), max(agent_id, id));
            }
        }
        if (from != to && path_table.table[to].size() >= t && path_table.table[from].size() > t) // edge conflicts
        {
            for (auto a1 : path_table.table[to][t - 1])
            {
                for (auto a2: path_table.table[from][t])
                {
                    if (a1 == a2)
                    {
                        succ = true;
                        colliding_pairs.emplace(min(agent_id, a1), max(agent_id, a1));
                        break;
                    }
                }
            }
        }
        //auto id = getAgentWithTarget(to, t);
        //if (id >= 0) // this agent traverses the target of another agent
        //    colliding_pairs.emplace(min(agent_id, id), max(agent_id, id));
        if (!path_table.goals.empty() && path_table.goals[to] < t) // target conflicts
        { // this agent traverses the target of another agent
            for (auto id : path_table.table[to][path_table.goals[to]]) // look at all agents at the goal time
            {
                if (agents[id].path.back().location == to) // if agent id's goal is to, then this is the agent we want
                {
                    succ = true;
                    colliding_pairs.emplace(min(agent_id, id), max(agent_id, id));
                    break;
                }
            }
        }
    }
    int goal = path.back().location; // target conflicts - some other agent traverses the target of this agent
    for (int t = (int)path.size(); t < path_table.table[goal].size(); t++)
    {
        for (auto id : path_table.table[goal][t])
        {
            succ = true;
            colliding_pairs.emplace(min(agent_id, id), max(agent_id, id));
        }
    }
    return succ;
}

void InitLNS::chooseDestroyHeuristicbyALNS()
{
    rouletteWheel();
    switch (selected_neighbor)
    {
        case 0 : init_destroy_strategy = TARGET_BASED; break;
        case 1 : init_destroy_strategy = COLLISION_BASED; break;
        case 2 : init_destroy_strategy = RANDOM_BASED; break;
        case 3 : init_destroy_strategy = SAT_BASED; break;
        default : cerr << "ERROR" << endl; exit(-1);
    }
}

bool InitLNS::generateNeighborByCollisionGraph()
{
    /*unordered_map<int, list<int>> G;
    for (int i = 0; i < (int)collision_graph.size(); i++)
    {
        if (!collision_graph[i].empty())
            G[i].assign(collision_graph[i].begin(), collision_graph[i].end());
    }
    assert(!G.empty());
    assert(neighbor_size <= (int)agents.size());
    set<int> neighbors_set;
    if ((int)G.size() < neighbor_size)
    {
        for (const auto& node : G)
            neighbors_set.insert(node.first);
        int count = 0;
        while ((int)neighbors_set.size() < neighbor_size && count < 10)
        {
            int a1 = *std::next(neighbors_set.begin(), rand() % neighbors_set.size());
            int a2 = randomWalk(a1);
            if (a2 != NO_AGENT)
                neighbors_set.insert(a2);
            else
                count++;
        }
    }
    else
    {
        int a = -1;
        while ((int)neighbors_set.size() < neighbor_size)
        {
            if (a == -1)
            {
                a = std::next(G.begin(), rand() % G.size())->first;
                neighbors_set.insert(a);
            }
            else
            {
                a = *std::next(G[a].begin(), rand() % G[a].size());
                auto ret = neighbors_set.insert(a);
                if (!ret.second) // no new element inserted
                    a = -1;
            }
        }
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by collision graph" << endl;
    return true;*/

    vector<int> all_vertices;
    all_vertices.reserve(collision_graph.size());
    for (int i = 0; i < (int)collision_graph.size(); i++)
    {
        if (!collision_graph[i].empty())
            all_vertices.push_back(i);
    }
    unordered_map<int, set<int>> G;
    auto v = all_vertices[rand() % all_vertices.size()]; // pick a random vertex
    findConnectedComponent(collision_graph, v, G);
    assert(G.size() > 1);

    assert(neighbor_size <= (int)agents.size());
    set<int> neighbors_set;
    if ((int)G.size() <= neighbor_size)
    {
        for (const auto& node : G)
            neighbors_set.insert(node.first);
        int count = 0;
        while ((int)neighbors_set.size() < neighbor_size && count < 10)
        {
            int a1 = *std::next(neighbors_set.begin(), rand() % neighbors_set.size());
            int a2 = randomWalk(a1);
            if (a2 != NO_AGENT)
                neighbors_set.insert(a2);
            else
                count++;
        }
    }
    else
    {
        int a = std::next(G.begin(), rand() % G.size())->first;
        neighbors_set.insert(a);
        while ((int)neighbors_set.size() < neighbor_size)
        {
            a = *std::next(G[a].begin(), rand() % G[a].size());
            neighbors_set.insert(a);
        }
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by collision graph" << endl;
    return true;

}
bool InitLNS::generateNeighborByTarget()
{
    int a = -1;
    auto r = rand() % (num_of_colliding_pairs * 2);
    int sum = 0;
    for (int i = 0 ; i < (int)collision_graph.size(); i++)
    {
        sum += (int)collision_graph[i].size();
        if (r <= sum and !collision_graph[i].empty())
        {
            a = i;
            break;
        }
    }
    assert(a != -1 and !collision_graph[a].empty());
    set<pair<int,int>> A_start; // an ordered set of (time, id) pair.
    set<int> A_target;


    for(int t = 0 ;t< path_table.table[agents[a].path_planner->start_location].size();t++){
        for(auto id : path_table.table[agents[a].path_planner->start_location][t]){
            if (id!=a)
                A_start.insert(make_pair(t,id));
        }
    }



    agents[a].path_planner->findMinimumSetofColldingTargets(goal_table,A_target);// generate non-wait path and collect A_target


    if (screen >= 3){
        cout<<"     Selected a : "<< a<<endl;
        cout<<"     Select A_start: ";
        for(auto e: A_start)
            cout<<"("<<e.first<<","<<e.second<<"), ";
        cout<<endl;
        cout<<"     Select A_target: ";
        for(auto e: A_target)
            cout<<e<<", ";
        cout<<endl;
    }

    set<int> neighbors_set;

    neighbors_set.insert(a);

    if(A_start.size() + A_target.size() >= neighbor_size-1){
        if (A_start.empty()){
            vector<int> shuffled_agents;
            shuffled_agents.assign(A_target.begin(),A_target.end());
            std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());
            neighbors_set.insert(shuffled_agents.begin(), shuffled_agents.begin() + neighbor_size-1);
        }
        else if (A_target.size() >= neighbor_size){
            vector<int> shuffled_agents;
            shuffled_agents.assign(A_target.begin(),A_target.end());
            std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());
            neighbors_set.insert(shuffled_agents.begin(), shuffled_agents.begin() + neighbor_size-2);

            neighbors_set.insert(A_start.begin()->second);
        }
        else{
            neighbors_set.insert(A_target.begin(), A_target.end());
            for(auto e : A_start){
                //A_start is ordered by time.
                if (neighbors_set.size()>= neighbor_size)
                    break;
                neighbors_set.insert(e.second);

            }
        }
    }
    else if (!A_start.empty() || !A_target.empty()){
        neighbors_set.insert(A_target.begin(), A_target.end());
        for(auto e : A_start){
            neighbors_set.insert(e.second);
        }

        set<int> tabu_set;
        while(neighbors_set.size()<neighbor_size){
            int rand_int = rand() % neighbors_set.size();
            auto it = neighbors_set.begin();
            std::advance(it, rand_int);
            a = *it;
            tabu_set.insert(a);

            if(tabu_set.size() == neighbors_set.size())
                break;

            vector<int> targets;
            for(auto p: agents[a].path){
                if(goal_table[p.location]>-1){
                    targets.push_back(goal_table[p.location]);
                }
            }

            if(targets.empty())
                continue;
            rand_int = rand() %targets.size();
            neighbors_set.insert(*(targets.begin()+rand_int));
        }
    }



    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by target" << endl;
    return true;
}
bool InitLNS::generateNeighborRandomly()
{
    if (neighbor_size >= agents.size())
    {
        neighbor.agents.resize(agents.size());
        for (int i = 0; i < (int)agents.size(); i++)
            neighbor.agents[i] = i;
        return true;
    }
    set<int> neighbors_set;
    auto total = num_of_colliding_pairs * 2 + agents.size();
    while(neighbors_set.size() < neighbor_size)
    {
        vector<int> r(neighbor_size - neighbors_set.size());
        for (auto i = 0; i < neighbor_size - neighbors_set.size(); i++)
            r[i] = rand() % total;
        std::sort(r.begin(), r.end());
        int sum = 0;
        for (int i = 0, j = 0; i < agents.size() and j < r.size(); i++)
        {
            sum += (int)collision_graph[i].size() + 1;
            if (sum >= r[j])
            {
                neighbors_set.insert(i);
                while (j < r.size() and sum >= r[j])
                    j++;
            }
        }
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors randomly" << endl;
    return true;
}

// Random walk; return the first agent that the agent collides with
int InitLNS::randomWalk(int agent_id)
{
    int t = rand() % agents[agent_id].path.size();
    int loc = agents[agent_id].path[t].location;
    while (t <= path_table.makespan and
           (path_table.table[loc].size() <= t or
           path_table.table[loc][t].empty() or
           (path_table.table[loc][t].size() == 1 and path_table.table[loc][t].front() == agent_id)))
    {
        auto next_locs = instance.getNeighbors(loc);
        next_locs.push_back(loc);
        int step = rand() % next_locs.size();
        auto it = next_locs.begin();
        loc = *std::next(next_locs.begin(), rand() % next_locs.size());
        t = t + 1;
    }
    if (t > path_table.makespan)
        return NO_AGENT;
    else
        return *std::next(path_table.table[loc][t].begin(), rand() % path_table.table[loc][t].size());
}

void InitLNS::writeIterStatsToFile(const string & file_name) const
{
    std::ofstream output;
    output.open(file_name);
    // header
    output << //"num of agents," <<
           "sum of costs," <<
           "num of colliding pairs," <<
           "runtime" << //"," <<
           //"MAPF algorithm" <<
           endl;

    for (const auto &data : iteration_stats)
    {
        output << //data.num_of_agents << "," <<
               data.sum_of_costs << "," <<
               data.num_of_colliding_pairs << "," <<
               data.runtime << //"," <<
               // data.algorithm <<
               endl;
    }
    output.close();
}

void InitLNS::writeResultToFile(const string & file_name, int sum_of_distances, double preprocessing_time) const
{
    std::ifstream infile(file_name);
    bool exist = infile.good();
    infile.close();
    if (!exist)
    {
        ofstream addHeads(file_name);
        addHeads << "runtime,num of collisions,solution cost,initial collisions,initial solution cost," <<
                 "sum of distances,iterations,group size," <<
                 "runtime of initial solution,area under curve," <<
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
            auc += prev->num_of_colliding_pairs * (curr->runtime - prev->runtime);
            prev = curr;
            ++curr;
        }
        auc += prev->num_of_colliding_pairs * (time_limit - prev->runtime);
    }

    ofstream stats(file_name, std::ios::app);
    stats << runtime << "," << iteration_stats.back().num_of_colliding_pairs << "," <<
          sum_of_costs << "," << iteration_stats.front().num_of_colliding_pairs << "," <<
          iteration_stats.front().sum_of_costs << "," << sum_of_distances << "," <<
          iteration_stats.size() << "," << average_group_size << "," <<
          iteration_stats.front().runtime << "," << auc << "," <<
          num_LL_expanded << "," << num_LL_generated << "," << num_LL_reopened << "," << num_LL_runs << "," <<
          preprocessing_time << "," << getSolverName() << "," << instance.getInstanceName() << endl;
    stats.close();
}

void InitLNS::printCollisionGraph() const
{
    cout << "Collision graph: ";
    int edges = 0;
    for (size_t i = 0; i < collision_graph.size(); i++)
    {
        for (int j : collision_graph[i])
        {
            if (i < j)
            {
                cout << "(" << i << "," << j << "),";
                edges++;
            }
        }
    }
    cout << endl <<  "|V|=" << collision_graph.size() << ", |E|=" << edges << endl;
}


unordered_map<int, set<int>>& InitLNS::findConnectedComponent(const vector<set<int>>& graph, int vertex,
                                                               unordered_map<int, set<int>>& sub_graph)
{
    std::queue<int> Q;
    Q.push(vertex);
    sub_graph.emplace(vertex, graph[vertex]);
    while (!Q.empty())
    {
        auto v = Q.front(); Q.pop();
        for (const auto & u : graph[v])
        {
            auto ret = sub_graph.emplace(u, graph[u]);
            if (ret.second) // insert successfully
                Q.push(u);
        }
    }
    return sub_graph;
}

void InitLNS::printPath() const
{
    for (const auto& agent : agents)
        cout << "Agent " << agent.id << ": " << agent.path << endl;
}

void InitLNS::printResult()
{
    average_group_size = - iteration_stats.front().num_of_agents;
    for (const auto& data : iteration_stats)
        average_group_size += data.num_of_agents;
    if (average_group_size > 0)
        average_group_size /= (double)(iteration_stats.size() - 1);
    assert(!iteration_stats.empty());
    cout << "\t" << getSolverName() << ": "
         << "runtime = " << runtime << ", "
         << "iterations = " << iteration_stats.size() << ", "
         << "colliding pairs = " << num_of_colliding_pairs << ", "
         << "initial colliding pairs = " << iteration_stats.front().num_of_colliding_pairs << ", "
         << "solution cost = " << sum_of_costs << ", "
         << "initial solution cost = " << iteration_stats.front().sum_of_costs << ", "
         << "failed iterations = " << num_of_failures << endl;
}

void InitLNS::clear()
{
    path_table.clear();
    collision_graph.clear();
    goal_table.clear();
}


bool InitLNS::validatePathTable() const
{
    for (auto i = 0; i < agents.size(); i++)
        assert(path_table.getPath(i) == &agents[i].path);
    return true;
}
