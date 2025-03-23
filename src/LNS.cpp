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
         num_of_iterations(num_of_iterations > 0 ? 0 : 2), // TODO: proč nefunguje?
         use_init_lns(use_init_lns),init_destory_name(init_destory_name),
         path_table(instance.map_size), pipp_option(pipp_option) {
    start_time = Time::now();
    replan_time_limit = time_limit / 100;
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
    else if (destory_name == "SAT") // new operator
        destroy_strategy = SAT;
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

// Vrací (sx, sy) odpovídající local_id v pořadí volných buněk v 'map' (2D pole s 1=volno, -1=prekážka)
pair<int, int> LNS::decodeLocalID(int local_id, const vector<vector<int>>& map) {
    int count = 0;
    int rows = map.size();
    int cols = map[0].size();
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (map[i][j] == 1) { // pokud je volná
                if (count == local_id)
                    return {i, j};
                count++;
            }
        }
    }
    // Pokud se index nepodaří nalézt, vrátíme chybový výsledek (můžete případně vyhodit výjimku)
    return {-1, -1};
}

/* getting the submap around one or more agents and identifying agents in these submaps */
pair<vector<vector<int>>, vector<int>> LNS::getSubmapAndAgents(int agent_id, int submap_size, int agent_location) {
    int map_width = 32;  // fixed for now...
    int map_height = 32; // fixed for now...

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
                    path_table.get_agents(conflicting_agents, global_pos); // collect all agents in the submap (at EVERY timestep)
                }
            }
        }
    }

    agents_in_submap.assign(conflicting_agents.begin(), conflicting_agents.end());
    return {submap, agents_in_submap};
}

void LNS::initializeSubmapData(const vector<vector<int>>& submap,
                               unordered_set<int>& submap_set,
                               unordered_map<int, pair<int, int>>& global_to_local) {
    cout << "Submap content (global positions):" << endl;
    for (size_t x = 0; x < submap.size(); ++x) {
        for (size_t y = 0; y < submap[x].size(); ++y) {
            int global_pos = submap[x][y];
            cout << global_pos << " ";
            if (global_pos != -1) {
                submap_set.insert(global_pos);
                global_to_local[global_pos] = {static_cast<int>(x), static_cast<int>(y)};
            }
        }
        cout << endl;
    }
}

vector<vector<int>> LNS::generateMapRepresentation(const vector<vector<int>>& submap,
                                                   const vector<int>& agents_in_submap,
                                                   int problematic_timestep) {
    cout << "Map content with agents and obstacles:" << endl;
    vector<vector<int>> map(submap.size(), vector<int>(submap[0].size(), 1));
    for (size_t x = 0; x < submap.size(); ++x) {
        for (size_t y = 0; y < submap[0].size(); ++y) {
            int global_pos = submap[x][y];
            if (global_pos == -1) cout << ". ";
            else if (instance.isObstacle(global_pos)) {
                cout << "X ";
                map[x][y] = -1;
            } else {
                bool is_agent = false;
                for (int agent : agents_in_submap) {
                    // TODO: agenti jsou zachytávání v problematic_timestep v submapě, ale ve skutečnosti tam mohou být v submapě v různé časy
                    // teoreticky nám to nevadí, ukazuje to místa, kde se agenti budou pohybovat v problematic_timestep
                    if (agents[agent].path[problematic_timestep].location == global_pos) {
                        cout << "A ";
                        is_agent = true;
                        break;
                    }
                }
                if (!is_agent) cout << ". ";
            }
        }
        cout << endl;
    }
    return map;
}

vector<int> LNS::getAgentsToReplan(const vector<int>& agents_in_submap,
                                   const unordered_set<int>& submap_set,
                                   int problematic_timestep)
{
    vector<int> agents_to_replan;
    cout << "\n[INFO] Identifikace agentů v submapě pro přeplánování "
         << "(kdy key_agent je tam v čase " << problematic_timestep << "):\n";

    for (int agent : agents_in_submap)
    {
        // TODO: co když je to ale druhý/třetí/... průchod agenta submapou, jehož interval už pokrývá problematic_timestep?
        // Najdeme PRVNÍ souvislý interval [t_min..t_max], ve kterém je agent v submapě.
        // Kdyby agent submapu opustil a pak se znovu vrátil, ignorujeme ten návrat.

        const auto& path = agents[agent].path;
        if (path.empty()) continue;

        int t_min = -1;
        int t_max = -1;

        // Projdeme celé path:
        //  1) hledáme první t, kdy je agent v submapě -> t_min
        //  2) dokud je agent v submapě, posouváme t_max
        //  3) jakmile agent submapu opustí, končíme (break)
        for (int t = 0; t < (int)path.size(); t++)
        {
            int loc = path[t].location;
            bool inSubmap = (submap_set.find(loc) != submap_set.end());

            if (t_min < 0)
            {
                // Ještě jsme žádný interval nezačali
                if (inSubmap)
                {
                    t_min = t;
                    t_max = t;
                }
                // jinak jen pokračujeme, dokud agent nevkročí do submapy
            }
            else
            {
                // Už jsme v intervalu
                if (inSubmap)
                    t_max = t; // agent stále uvnitř submapy
                else break; // Agent poprvé opustil submapu => skončíme
            }
        }

        // Pokud t_min zůstalo -1, agent do submapy vůbec nevkročil
        if (t_min == -1) {
            // vynecháváme
            continue;
        }

        // Tady máme interval [t_min..t_max], kde agent poprvé pobývá v submapě
        // a nepouštíme se do dalších potenciálních návratů

        // Ověříme, zda [t_min..t_max] obsahuje problematický čas
        if (t_min <= problematic_timestep && problematic_timestep <= t_max)
        {
            agents_to_replan.push_back(agent);
            cout << "  Agent " << agent
                 << " (první interval v submapě je ["
                 << t_min << ".." << t_max << "]),"
                 << " pokrývá i čas " << problematic_timestep
                 << ", přidán k přeplánování.\n";
        }
    }

    if (agents_to_replan.empty())
        cout << "[INFO] Žádný agent nebyl v submapě ve stejnou chvíli (čas "
             << problematic_timestep << ").\n";

    return agents_to_replan;
}

// Vrátí mapu [agent -> vektor dvojic (sx, sy) v submapě],
// od T_sync do doby, kdy agent submapu opustí.
unordered_map<int, vector<pair<int,int>>>
LNS::findLocalPaths(const vector<int>& agents_to_replan,
                    const vector<vector<int>>& submap,
                    const unordered_set<int>& submap_set,
                    const unordered_map<int, pair<int,int>>& global_to_local,
                    int T_sync) {
    // TODO: co když agent prochází submapou opakovaně? Musíme vybrat ten souvislý úsek lokální cesty, která má v intervalu problematic_timestep
    // pro každého agent vytvoříme sekvenci (sx, sy) lokálních souřadnic
    // pokud agent není v submapě v čase T_sync, vynecháme ho
    // jakmile agent submapu opustí, končíme

    unordered_map<int, vector<pair<int,int>>> local_paths;

    cout << "\n[INFO] Tvorba lokálních cest (sx, sy) v submapě pro T_sync = "
         << T_sync << endl;

    for (int agent : agents_to_replan)
    {
        // 1) Zajistíme, že agent má definovanou pozici v T_sync
        if ((size_t)T_sync >= agents[agent].path.size()) {
            cout << "[WARN] Agent " << agent
                 << " nemá definovanou pozici v čase T_sync="
                 << T_sync << ". Přeskakuji.\n";
            continue;
        }

        // 2) Ověříme, že je agent v submapě v T_sync
        int loc_at_Tsync = agents[agent].path[T_sync].location;
        if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
            cout << "[WARN] Agent " << agent
                 << " není v submapě v čase " << T_sync
                 << ". Přeskakuji.\n";
            continue;
        }

        // 3) Najdeme poslední čas, dokdy agent v submapě zůstává
        int last_time_in_submap = -1;
        for (int t = T_sync; t < (int)agents[agent].path.size(); t++) {
            int glob_loc = agents[agent].path[t].location;
            // Dokud je glob_loc v submapě, posouváme last_time_in_submap
            if (submap_set.find(glob_loc) != submap_set.end())
                last_time_in_submap = t;
            else break; // agent submapu opustil
        }

        if (last_time_in_submap == -1) {
            // Teoreticky by se to nemělo stát,
            // protože loc_at_Tsync je v submapě
            cout << "[WARN] Agent " << agent
                 << " v submapě vlastně není? (podivné)\n";
            continue;
        }

        // 4) Postavíme reálnou lokální cestu v (sx, sy)
        //    od T_sync do last_time_in_submap
        vector<pair<int,int>> path_local;
        for (int t = T_sync; t <= last_time_in_submap; t++)
        {
            int glob_loc = agents[agent].path[t].location;
            auto it = global_to_local.find(glob_loc);
            if (it == global_to_local.end()) {
                // Mělo by se stávat jen pokud agent
                // reálně vyběhl z submapy
                cout << "[ERROR] Agent " << agent
                     << " v case " << t
                     << " je mimo submapu, ale last_time_in_submap="
                     << last_time_in_submap << endl;
                break;
            }
            // (sx, sy) = lokální souřadnice v submapě
            int sx = it->second.first;
            int sy = it->second.second;

            path_local.emplace_back(sx, sy);
        }

        // 5) uložíme do mapy
        local_paths[agent] = path_local;

        cout << "  Agent " << agent
             << " (globální cesty od T=" << T_sync
             << " do " << last_time_in_submap << ") má lokální dráhu: ";

        for (auto& [sx, sy] : path_local)
            cout << "(" << sx << "," << sy << ") ";
        cout << endl;
    }

    return local_paths;
}

bool LNS::solveWithSAT(
        vector<vector<int>>& map,
        const unordered_map<int, vector<pair<int,int>>>& local_paths,
        vector<int>& agents_to_replan,
        const vector<vector<int>>& submap,
        int T_sync) {
    cout << "\n[DEBUG] Kontrola vstupních dat pro SAT solver:\n";
    cout << "  - Počet agentů k přeplánování: " << agents_to_replan.size() << endl;

    // 1) Sestavení start/goal pro solver (z local_paths)
    vector<pair<int,int>> start_positions;
    vector<pair<int,int>> goal_positions;
    std::map<int,int> original_local_lengths;

    for (int agent : agents_to_replan) {
        auto it = local_paths.find(agent);
        if (it == local_paths.end() || it->second.empty()) {
            cout << "[WARN] Agent " << agent
                 << " nemá local_path => vynecháváme.\n";
            continue;
        }
        const auto& path = it->second; // vector<pair<int,int>>
        // start = path.front(), goal = path.back()
        start_positions.push_back(path.front());
        goal_positions.push_back(path.back());
        original_local_lengths[agent] = (int) path.size();

        cout << "[DEBUG] Agent " << agent
             << " má původní lokální dráhu délky: "
             << path.size()
             << " => Start ("
             << path.front().first << "," << path.front().second << "), "
             << "Goal ("
             << path.back().first << "," << path.back().second << ")\n";
    }

    // 2) Vytvoření SAT instance
    auto inst = std::make_unique<_MAPFSAT_Instance>(map, start_positions, goal_positions);
    auto solver = std::make_unique<_MAPFSAT_DisappearAtGoal>();
    auto log    = std::make_unique<_MAPFSAT_Logger>(inst.get(), "disappear_at_goal", 2);

    cout << "SAT instance and solver created.\n";

    solver->SetData(inst.get(), log.get(), 300, "", false, true);
    inst->SetAgents((int)start_positions.size());
    log->NewInstance((int)start_positions.size());

    // 3) Spustíme solver
    int result = solver->Solve((int)start_positions.size(), 0, true);
    cout << "Solver returned: " << result << endl;

    if (result != 0) {
        cout << "SAT solver failed.\n";
        return false;
    }

    // 4) Získáme nové cesty od solveru (v indexech submapy)
    // TODO: odseknout opakující se konečné pozice
    vector<vector<int>> plan = solver->GetPlan();

    for (auto& path_for_agent : plan) {
        if (path_for_agent.empty()) continue;
        // Smaž opakované indexy na úplném konci (dokud se opakují v cíli):
        while (path_for_agent.size() > 1 && path_for_agent.back() == path_for_agent[path_for_agent.size() - 2]) {
            path_for_agent.pop_back();
        }
    }

    // Debug: vypsat novou lokální cestu (v 1D indexech submapy).
    for (size_t a = 0; a < plan.size(); ++a) {
        cout << "[DEBUG] Agent (index) " << agents_to_replan[a]
             << " | Nová lokální cesta (submap idx): ";
        for (auto lid : plan[a]) {
            cout << lid << " ";
        }
        cout << endl;
    }

    // 5) Aktualizujeme updated_path pro každého agenta (prefix + nová lokální dráha + sufix)
    for (size_t a = 0; a < plan.size(); ++a) {
        int agent_id = agents_to_replan[a];

        // (A) Původní délka lokální části
        int old_local_length = original_local_lengths[agent_id];
        // (B) Nová délka lokální části
        int new_local_length = (int) plan[a].size();

        cout << "[INFO] Aktualizace cesty pro agenta " << agent_id
             << " | Původní lokální délka: " << old_local_length
             << " | Nová lokální délka: " << new_local_length << endl;

        // (1) Zkopírujeme prefix do T_sync
        vector<PathEntry> updated_path(
            agents[agent_id].path.begin(),
            agents[agent_id].path.begin() + T_sync
        );

        // (2) Vložíme novou lokální trasu (dekódovanou do glob. ID)
        for (int t = 0; t < new_local_length; t++)
        {
            int local_id = plan[a][t];
            // decodeLocalID => (sx,sy) submap
            // *** POZOR na správné x/y a submap[x][y] vs submap[y][x]! ***
            pair<int, int> coords = decodeLocalID(local_id, map);
            int sx = coords.first;
            int sy = coords.second;
            if (sx == -1 || sy == -1) {
                cout << "[ERROR] agent " << agent_id
                     << " local_id=" << local_id
                     << " nelze dekódovat do platných souřadnic.\n";
                continue;
            }
            // Globální ID z submapy:
            int global_id = submap[sx][sy];
            cout << "[DEBUG] agent " << agent_id << " t=" << t
                 << " => decoded (sx,sy)=(" << sx << "," << sy
                 << ") => global_id=" << global_id << endl;
            updated_path.push_back(PathEntry(global_id));
        }

        // (3) Připojíme suffix. Původní suffix začínal na indexu T_sync + old_local_length
        // (tj. tam končila stará lokální část). V nové cestě jsme skončili na indexu T_sync + new_local_length - 1
        // => suffix nalepíme od původního suffix_start dále.
        int old_suffix_start = T_sync + old_local_length;
        if (old_suffix_start < (int)agents[agent_id].path.size())
        {
            for (int t = old_suffix_start; t < (int)agents[agent_id].path.size(); t++) {
                updated_path.push_back(PathEntry(agents[agent_id].path[t].location));
            }
        }

        // (4) Volitelná kontrola navázání prefix->lokální
        if (T_sync > 0 && (size_t)T_sync < updated_path.size()) {
            int prefix_last = agents[agent_id].path[T_sync].location; // TODO: nevím, jestli dává smysl kontrola návaznosti prefixu
            int local_first = updated_path[T_sync].location;
            if (prefix_last != local_first) {
                cout << "[WARN] agent " << agent_id
                     << " prefix->lokální navaznost se liší: "
                     << prefix_last << " != " << local_first << endl;
            }
        }

        cout << "[INFO] Původní délka cesty agenta " << agent_id
             << " je: " <<  agents[agent_id].path.size() << endl;
        // (5) Uložíme hotovou cestu
        agents[agent_id].path = updated_path;
        cout << "[INFO] Cesta pro agenta " << agent_id
             << " aktualizována, výsledná délka: "
             << updated_path.size() << endl;
    }

    cout << "Paths successfully updated.\n";
    return true;
}

/* NOTE:
 * náš operátor by se měl správně volat na vyřešení konfliktů, jakmile se konflikty vyřeší,
 * bude se volat optimalizace (na to máme nějaký určený čas, třeba 10 s), nicméně náš operátor při optimalizaci
 * může působit nové konflikty, což ostatní operátory nedělají, proto je potřeba najít v projektu flag,
 * který přepíná z řešení konfliktů na optimalizaci, u našeho operátoru by bylo třeba přepínat mezi
 * naším optimalizačním operátorem a řešením konfliktů (na přeskáčku přepínat: je konflikt -> vyřeš -> optimalizace -> ...)
 * */
// --------------------------------------------------------
// DESTROY fáze: generateNeighborBySAT() – najde submapu, agenty, T_sync atd.
// --------------------------------------------------------
bool LNS::generateNeighborBySAT() {
    cout << "====================" << endl;
    cout << "SAT destroy operator called." << endl;

    auto [key_agent_id, problematic_timestep] = findMostDelayedAgent();
    if (key_agent_id < 0) {
        cout << "No delayed agent found." << endl;
        return false;
    }

    int agent_loc = agents[key_agent_id].path[problematic_timestep].location; // globalID of the cell in 1D matrix
    int submap_size = 25;
    int map_width = 32;  // případně dynamicky z instance
    int map_height = 32; // případně dynamicky z instance

    // int agent_id, int submap_size, int agent_location
    auto [submap, agents_in_submap] = getSubmapAndAgents(key_agent_id, submap_size, agent_loc);

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
    ignored_agents.insert(key_agent_id);

    return true;
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
        cout << "[WARN] SAT solver failed to find a valid solution." << endl;
        return false;
    }

    // úspěch – spočítáme novou sum_of_costs
    neighbor.sum_of_costs = 0;
    for (int ag : agents_to_replan) {
        neighbor.sum_of_costs += (int)agents[ag].path.size() - 1;
    }

    if (neighbor.sum_of_costs <= neighbor.old_sum_of_costs) {
        // akceptujeme novou cestu
        for (int ag : agents_to_replan) {
            path_table.insertPath(ag, agents[ag].path);
        }
        return true;
    } else {
        // revert
        cout << "[INFO] New SAT solution is worse, reverting." << endl;
        for (int i = 0; i < (int)neighbor.agents.size(); i++) {
            int a = neighbor.agents[i];
            path_table.deletePath(agents[a].id, agents[a].path);
            agents[a].path = neighbor.old_paths[i];
            path_table.insertPath(agents[a].id, agents[a].path);
        }
        neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        return false;
    }
}

bool LNS::run()
{
    // only for statistic analysis, and thus is not included in runtime
    sum_of_distances = 0;
    for (const auto & agent : agents)
    {
        sum_of_distances += agent.path_planner->my_heuristic[agent.path_planner->start_location];
    }

    initial_solution_runtime = 0;
    start_time = Time::now();
    bool succ = getInitialSolution();
    initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
    if (!succ && initial_solution_runtime < time_limit)
    {
        if (use_init_lns)
        {
            init_lns = new InitLNS(instance, agents, time_limit - initial_solution_runtime,
                    replan_algo_name,init_destory_name, neighbor_size, screen);
            succ = init_lns->run();
            if (succ) // accept new paths
            {
                path_table.reset();
                for (const auto & agent : agents)
                {
                    path_table.insertPath(agent.id, agent.path);
                }
                init_lns->clear();
                initial_sum_of_costs = init_lns->sum_of_costs;
                sum_of_costs = initial_sum_of_costs;
            }
            initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
        }
        else // use random restart
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
    if (succ)
    {
        if (screen >= 1)
            cout << "Initial solution cost = " << initial_sum_of_costs << ", "
                 << "runtime = " << initial_solution_runtime << endl;
    }
    else
    {
        cout << "Failed to find an initial solution in "
             << runtime << " seconds after  " << restart_times << " restarts" << endl;
        return false; // terminate because no initial solution is found
    }

    while (runtime < time_limit && iteration_stats.size() <= num_of_iterations)
    {
        cout.flush();
        runtime =((fsec)(Time::now() - start_time)).count();
        if(screen >= 1)
            validateSolution(); // vrací exit(-1) při nalezení konfliktu -> problém
        if (ALNS)
            chooseDestroyHeuristicbyALNS();

        switch (destroy_strategy)
        {
            case RANDOMWALK:
                succ = generateNeighborByRandomWalk();
                break;
            case INTERSECTION:
                succ = generateNeighborByIntersection();
                break;
            case RANDOMAGENTS:
                neighbor.agents.resize(agents.size());
                for (int i = 0; i < (int)agents.size(); i++)
                    neighbor.agents[i] = i;
                if (neighbor.agents.size() > neighbor_size)
                {
                    std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
                    neighbor.agents.resize(neighbor_size);
                }
                succ = true;
                break;
            case SAT: { // destroy i repair dohromady
                // Jedna iterace => až MAX_SAT_ATTEMPTS pro SAT
                const int MAX_SAT_ATTEMPTS = 10;
                bool sat_success = false;
                for (int attempt = 0; attempt < MAX_SAT_ATTEMPTS && !sat_success; attempt++) {
                    if (!generateNeighborBySAT())
                        continue; // nepodařilo se najít validní neighborhood
                    // Máme neighborhood, teď se pokusíme přeplánovat pomocí SAT
                    if (runSAT())
                        sat_success = true; // SAT operátor našel validní řešení
                }
                succ = sat_success;
            } break;
            default:
                cerr << "Wrong neighbor generation strategy" << endl;
                exit(-1);
        }

        if(!succ)
            continue;

        // store the neighbor information
        neighbor.old_paths.resize(neighbor.agents.size());
        neighbor.old_sum_of_costs = 0;
        for (int i = 0; i < (int)neighbor.agents.size(); i++)
        {
            if (replan_algo_name == "PP")
                neighbor.old_paths[i] = agents[neighbor.agents[i]].path;
            path_table.deletePath(neighbor.agents[i], agents[neighbor.agents[i]].path);
            neighbor.old_sum_of_costs += agents[neighbor.agents[i]].path.size() - 1;
        }

        if (replan_algo_name == "EECBS")
            succ = runEECBS();
        else if (replan_algo_name == "CBS")
            succ = runCBS();
        else if (replan_algo_name == "PP")
            succ = runPP();
        else if (replan_algo_name == "SAT") { /* nic */ }
        else
        {
            cerr << "Wrong replanning strategy" << endl;
            exit(-1);
        }

        if (ALNS) // update destroy heuristics
        {
            if (neighbor.old_sum_of_costs > neighbor.sum_of_costs )
                destroy_weights[selected_neighbor] =
                        reaction_factor * (neighbor.old_sum_of_costs - neighbor.sum_of_costs) / neighbor.agents.size()
                        + (1 - reaction_factor) * destroy_weights[selected_neighbor];
            else
                destroy_weights[selected_neighbor] =
                        (1 - decay_factor) * destroy_weights[selected_neighbor];
        }
        runtime = ((fsec)(Time::now() - start_time)).count();
        sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;
        if (screen >= 1)
            cout << "Iteration " << iteration_stats.size() << ", "
                 << "group size = " << neighbor.agents.size() << ", "
                 << "solution cost = " << sum_of_costs << ", "
                 << "remaining time = " << time_limit - runtime << endl;
        iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, replan_algo_name);
    }


    average_group_size = - iteration_stats.front().num_of_agents;
    for (const auto& data : iteration_stats)
        average_group_size += data.num_of_agents;
    if (average_group_size > 0)
        average_group_size /= (double)(iteration_stats.size() - 1);

    cout << getSolverName() << ": "
         << "runtime = " << runtime << ", "
         << "iterations = " << iteration_stats.size() << ", "
         << "solution cost = " << sum_of_costs << ", "
         << "initial solution cost = " << initial_sum_of_costs << ", "
         << "failed iterations = " << num_of_failures << endl;
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
        succ = runSAT();
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

pair<int, int> LNS::findMostDelayedAgent() {
    int max_delays = -1;
    int agent_with_max_delays = -1;
    int most_problematic_timestep = -1;

    for (const auto& agent : agents) {
        // přeskoč agenty, kteří už selhali
        // TODO: tohle není nejmoudřejší, protože je nechceme přeskakovat natrvalo, ale prozatím to nevadí
        // TODO: nikdy nevymazáváme ignored_agents.clear()
        if (ignored_agents.find(agent.id) != ignored_agents.end())
            continue;

        auto [agent_max_delays, problematic_timestep] = agent.getMostProblematicDelay(path_table);

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
}

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
