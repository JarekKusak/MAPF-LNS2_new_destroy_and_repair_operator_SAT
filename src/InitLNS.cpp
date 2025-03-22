#include "InitLNS.h"
#include <queue>
#include <algorithm>
#include "GCBS.h"
#include "PBS.h"

#include "../include/MAPF.hpp"

InitLNS::InitLNS(const Instance& instance, vector<Agent>& agents, double time_limit,
         const string & replan_algo_name, const string & init_destory_name, int neighbor_size, int screen) :
         BasicLNS(instance, time_limit, neighbor_size, screen), agents(agents), replan_algo_name(replan_algo_name),
         path_table(instance.map_size, agents.size()), collision_graph(agents.size()), goal_table(instance.map_size, -1)
 {
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
     else
     {
         cerr << "Init Destroy heuristic " << init_destory_name << " does not exists. " << endl;
         exit(-1);
     }

     for (auto& i:agents) {
         goal_table[i.path_planner->goal_location] = i.id;
     }

 }

/**
* findConflictAgent
*  - Najde agenta s největším počtem konfliktů v collision_graph (tj. největší size() v collision_graph[a]).
*  - Pokud collision_graph[a].size() == 0 pro všechny a, vrátí {-1,-1} => žádné konflikty.
*  - Pro zjednodušení vracíme conflict_time=0;
*    v reálné implementaci ho lze najít analýzou path_tableWC (např. hasCollisions / getLastCollisionTimestep).
*/
pair<int,int> InitLNS::findConflictAgent()
{
    int agent_with_most_conflicts = -1;
    size_t max_conflicts = 0;

    // Pokud collision_graph.size() < agents.size(),
    //   pak je tam buď chybějící agent, nebo zatím nebyla data aktualizována.
    //   Zde předpokládáme, že je vše sladěné.
    for (int a = 0; a < (int)collision_graph.size(); a++)
    {
        // collision_graph[a] je množina agentů, se kterými je agent a v konfliktu
        size_t conflicts = collision_graph[a].size();
        if (conflicts > max_conflicts)
        {
            max_conflicts = conflicts;
            agent_with_most_conflicts = a;
        }
    }

    // Pokud agent_with_most_conflicts zůstane -1 nebo max_conflicts == 0 => žádné konflikty
    if (agent_with_most_conflicts < 0 || max_conflicts == 0)
    {
        return {-1, -1}; // signál, že nejsou konflikty
    }

    // Tady byste mohli (volitelně) najít konkrétní čas konfliktu
    //   pro agent_with_most_conflicts, buď analýzou path_tableWC
    //   nebo jinými mechanismy. Zatím vracíme 0 jako placeholder.
    int conflict_time = 0;

    return { agent_with_most_conflicts, conflict_time };
}


// Vrací (sx, sy) odpovídající local_id v pořadí volných buněk v 'map' (2D pole s 1=volno, -1=prekážka)
pair<int, int> InitLNS::decodeLocalID(int local_id, const vector<vector<int>>& map) {
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

/* getting the submap around one or more agents and identifying agents in these submaps
pair<vector<vector<int>>, vector<int>> InitLNS::getSubmapAndAgents(int agent_id, int submap_size, int agent_location) {
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

void InitLNS::initializeSubmapData(const vector<vector<int>>& submap,
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

vector<vector<int>> InitLNS::generateMapRepresentation(const vector<vector<int>>& submap,
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

vector<int> InitLNS::getAgentsToReplan(const vector<int>& agents_in_submap,
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
InitLNS::findLocalPaths(const vector<int>& agents_to_replan,
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

bool InitLNS::solveWithSAT(
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
 *
// --------------------------------------------------------
// DESTROY fáze: generateNeighborBySAT() – najde submapu, agenty, T_sync atd.
// --------------------------------------------------------
bool InitLNS::generateNeighborBySAT() {
    cout << "====================" << endl;
    cout << "SAT operator called." << endl;

    //recently_replanned_agents.clear();

    const int MAX_ATTEMPTS = 10;

    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++)
    {
        auto [key_agent_id, problematic_timestep] = findMostDelayedAgent();
        if (key_agent_id < 0) {
            cout << "No delayed agent found." << endl;
            return false; // return true?
        }

        int agent_loc = agents[key_agent_id].path[problematic_timestep].location; // globalID of the cell in 1D matrix
        int submap_size = 25;
        auto [submap, agents_in_submap] = getSubmapAndAgents(key_agent_id, submap_size, agent_loc);

        unordered_set<int> submap_set;
        unordered_map<int, pair<int, int>> global_to_local;
        initializeSubmapData(submap, submap_set, global_to_local); // +submap -submap_set -global_to_local

        vector<vector<int>> map = generateMapRepresentation(submap, agents_in_submap, problematic_timestep);

        // NOTE
        // sum of cost známe, spočítáme optimální sum of cost od startu do cíle (DFS) každého agenta
        // podle toho můžeme říkat solveru, aby našel řešení s nějakou danou cenou (inkrementálně navyšujem) - delta
        // pokud to bude vyšší než to co známe, tak zahodíme submapu (nemá cenu řešit)
        // avoid existuje

        // TODO: avoid na agenty, kteří submapou prochází v budoucnosti pomocí path_table.getAgents
        vector<int> agents_to_replan = getAgentsToReplan(agents_in_submap, submap_set, problematic_timestep);
        if (agents_to_replan.empty()) return false; // return true?

        int T_sync = problematic_timestep; // bude se synchronizovat podle nejproblematičtějšího agenta

        // =================== DEBUG ======================
        // =================== DEBUG ======================
        // =================== DEBUG ======================
        vector<pair<int,int>> start_positions, goal_positions;
        for (int agent : agents_to_replan)
        {
            int start_global = -1, goal_global = -1;
            int start_time = T_sync, goal_time = -1;  // start je T_sync
            // 1) Ověříme, že agent je definován v T_sync
            if ((size_t)T_sync >= agents[agent].path.size()) {
                cout << "[ERROR] Agent " << agent
                     << " nemá definovanou pozici v čase T_sync!\n";
                continue;
            }

            // 2) Zjistíme, zda je v submapě v T_sync
            int loc_at_Tsync = agents[agent].path[T_sync].location;
            if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
                cout << "[WARNING] Agent " << agent
                     << " není v submapě v čase T_sync!\n";
                continue;
            }
            start_global = loc_at_Tsync;

            // 3) Najdeme první okamžik, kdy agent submapu opouští.
            //    => goal_time bude poslední t, pro který agent byl uvnitř submapy
            goal_time = -1;
            for (int t = T_sync; t < (int)agents[agent].path.size(); t++) {
                int location = agents[agent].path[t].location;
                if (submap_set.find(location) != submap_set.end()) {
                    // agent je pořád v submapě
                    goal_global = location;
                    goal_time = t;
                } else break; // agent submapu opustil poprvé => končíme
            }

            if (goal_time == -1) {
                // pokud by se stalo, že agent nebyl v submapě ani v T_sync,
                // ale to se nestane, protože výše jsme kontrolovali loc_atTsync
                cout << "[WARNING] Agent " << agent
                     << " nemá platnou cílovou pozici (už T_sync je mimo?)!\n";
                continue;
            }

            // 4) Převod globálních souřadnic start/goal na lokální (sx, sy)
            auto itS = global_to_local.find(start_global);
            auto itG = global_to_local.find(goal_global);

            if (itS == global_to_local.end()) {
                cout << "[ERROR] Startovní pozice agenta " << agent
                     << " není v global_to_local!\n";
                continue;
            }
            if (itG == global_to_local.end()) {
                cout << "[ERROR] Cílová pozice agenta " << agent
                     << " není v global_to_local!\n";
                continue;
            }

            // Uložit do start_positions, goal_positions pro debug
            start_positions.push_back(itS->second);
            goal_positions.push_back(itG->second);

            // 5) Výpis
            cout << "Agent " << agent
                 << " | Start (globální): " << start_global
                 << " → (lokální): (" << itS->second.first << ", " << itS->second.second << ")"
                 << " v čase " << T_sync
                 << " | První opuštění submapy v t=" << (goal_time+1)
                 << " => Cíl (globální): " << goal_global
                 << " → (lokální): (" << itG->second.first << ", " << itG->second.second << ")"
                 << " v čase " << goal_time << endl;
        }
        // =================== DEBUG ======================
        // =================== DEBUG ======================
        // =================== DEBUG ======================

        // Namísto volání findLocalPaths + solveWithSAT přímo tady
        // jen uložíme data do neighbor.* pro pozdější "repair" fází:
        neighbor.agents = agents_to_replan;
        neighbor.submap = submap;      // Uložíme si pro runSAT()
        neighbor.submap_set = submap_set;
        neighbor.global_to_local = global_to_local;
        neighbor.map = map;
        neighbor.T_sync = T_sync;

        // Také si můžeme poznamenat, že jsme si "vybrali" tento submap pro replan
        // => vrátíme true, a v run() potom dojde k volbě "runSAT()" namísto runPP apod.
        // (případně, pokud chceme vícekrokové pokusy, ponecháme cycle).
        return true; // zničení se povedlo, vrátíme se, abychom v run() spustili "repair" (runSAT).

        // Pokud by to selhalo, pak:
        // ignored_agents.insert(key_agent_id);
        // atd. – ale pro ukázku to teď vypneme.
    }

    cout << "[ERROR] SAT se nepodařilo aplikovat na žádného agenta po "
         << MAX_ATTEMPTS << " pokusech.\n";
    return false;
}

// --------------------------------------------------------
// REPAIR fáze: runSAT() – zavolá findLocalPaths + solveWithSAT,
//              a upraví cesty agentů + path_table
// --------------------------------------------------------
bool InitLNS::runSAT()
{
    cout << "====================" << endl;
    cout << "[REPAIR] SAT operator – spouštím subproblém." << endl;

    // Získáme z neighbor.* všechny uložené informace
    const auto& agents_to_replan = neighbor.agents;
    const auto& submap           = neighbor.submap;
    const auto& submap_set       = neighbor.submap_set;
    const auto& global_to_local  = neighbor.global_to_local;
    const auto& map              = neighbor.map;
    int T_sync                   = neighbor.T_sync;

    // Zavoláme původní findLocalPaths
    auto local_paths = findLocalPaths(agents_to_replan,
                                      submap,
                                      submap_set,
                                      global_to_local,
                                      T_sync);

    // Nyní spustíme solveWithSAT (původní kód).
    bool success = solveWithSAT(
            const_cast<vector<vector<int>>&>(map), // solveWithSAT bere non-const
            local_paths,
            const_cast<vector<int>&>(agents_to_replan),
            submap,
            T_sync);

    if (!success) {
        cout << "[WARN] SAT solver failed to find a valid solution." << endl;
        // Vrátíme false => run() to pozná a zahodí neighbor
        return false;
    }

    // Úspěch – spočítáme novou sum_of_costs
    neighbor.sum_of_costs = 0;
    for (int ag : agents_to_replan) {
        neighbor.sum_of_costs += (int)agents[ag].path.size() - 1;
    }

    // Porovnáme s old_sum_of_costs => buď akceptujeme, nebo revert
    if (neighbor.sum_of_costs <= neighbor.old_sum_of_costs) {
        // Akceptujeme novou cestu => zaneseme do path_table
        for (int ag : agents_to_replan) {
            path_table.insertPath(ag, agents[ag].path);
        }
        return true;
    }
    else {
        cout << "[INFO] New SAT solution is worse, reverting." << endl;
        // Vrátíme agentům staré cesty
        for (int i = 0; i < (int)neighbor.agents.size(); i++) {
            int a = neighbor.agents[i];
            path_table.deletePath(agents[a].id, agents[a].path);
            agents[a].path = neighbor.old_paths[i];
            path_table.insertPath(agents[a].id, agents[a].path);
        }
        neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        return false;
    }
}*/

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
    if (runtime >= time_limit and !succ)
    {
        printResult();
        return false;
    }

    vector<Path*> paths(agents.size());
    for (auto i = 0; i < agents.size(); i++)
        paths[i] = &agents[i].path;
    while (runtime < time_limit and num_of_colliding_pairs > 0)
    {
        assert(instance.validateSolution(paths, sum_of_costs, num_of_colliding_pairs));
        if (ALNS)
            chooseDestroyHeuristicbyALNS();

        switch (init_destroy_strategy)
        {
            case TARGET_BASED:
                succ = generateNeighborByTarget();
                break;
            case COLLISION_BASED:
                succ = generateNeighborByCollisionGraph();
                break;
            case RANDOM_BASED:
                succ = generateNeighborRandomly();
                break;
            default:
                cerr << "Wrong neighbor generation strategy" << endl;
                exit(-1);
        }
        if(!succ || neighbor.agents.empty())
            continue;

        // get colliding pairs
        neighbor.old_colliding_pairs.clear();
        for (int a : neighbor.agents)
        {
            for (auto j: collision_graph[a])
            {
                neighbor.old_colliding_pairs.emplace(min(a, j), max(a, j));
            }
        }
        if (neighbor.old_colliding_pairs.empty()) // no need to replan
        {
            assert(init_destroy_strategy == RANDOM_BASED);
            if (ALNS) // update destroy heuristics
            {
                destroy_weights[selected_neighbor] = (1 - decay_factor) * destroy_weights[selected_neighbor];
            }
            continue;
        }

        // store the neighbor information
        neighbor.old_paths.resize(neighbor.agents.size());
        neighbor.old_sum_of_costs = 0;
        for (int i = 0; i < (int)neighbor.agents.size(); i++)
        {
            int a = neighbor.agents[i];
            if (replan_algo_name == "PP" || neighbor.agents.size() == 1)
                neighbor.old_paths[i] = agents[a].path;
            path_table.deletePath(neighbor.agents[i]);
            neighbor.old_sum_of_costs += (int) agents[a].path.size() - 1;
        }
        if (screen >= 2)
        {
            cout << "Neighbors: ";
            for (auto a : neighbor.agents)
                cout << a << ", ";
            cout << endl;
            cout << "Old colliding pairs (" << neighbor.old_colliding_pairs.size() << "): ";
            for (const auto & p : neighbor.old_colliding_pairs)
            {
                cout << "(" << p.first << "," << p.second << "), ";
            }
            cout << endl;

        }

        if (replan_algo_name == "PP" || neighbor.agents.size() == 1)
            succ = runPP();
        else if (replan_algo_name == "GCBS")
            succ = runGCBS();
        else if (replan_algo_name == "SAT") // TODO: dočasně
            succ = runPP();
        else if (replan_algo_name == "PBS")
            succ = runPBS();
        else
        {
            cerr << "Wrong replanning strategy" << endl;
            exit(-1);
        }

        if (ALNS) // update destroy heuristics
        {
            if (neighbor.colliding_pairs.size() < neighbor.old_colliding_pairs.size())
                destroy_weights[selected_neighbor] =
                        reaction_factor * (double)(neighbor.old_colliding_pairs.size() -
                        neighbor.colliding_pairs.size()) // / neighbor.agents.size()
                        + (1 - reaction_factor) * destroy_weights[selected_neighbor];
            else
                destroy_weights[selected_neighbor] =
                        (1 - decay_factor) * destroy_weights[selected_neighbor];
        }
        if (screen >= 2)
            cout << "New colliding pairs = " << neighbor.colliding_pairs.size() << endl;
        if (succ) // update collision graph
        {
            num_of_colliding_pairs += (int)neighbor.colliding_pairs.size() - (int)neighbor.old_colliding_pairs.size();
            for(const auto& agent_pair : neighbor.old_colliding_pairs)
            {
                collision_graph[agent_pair.first].erase(agent_pair.second);
                collision_graph[agent_pair.second].erase(agent_pair.first);
            }
            for(const auto& agent_pair : neighbor.colliding_pairs)
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
        int id = *p;
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
