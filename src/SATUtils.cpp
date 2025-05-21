#include "SATUtils.h"

#include <iostream>
#include <algorithm>

namespace SATUtils {
    std::pair<int, int> decodeLocalID(int local_id, const std::vector<std::vector<int>>& map) {
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
        // Pokud se index nepodaří nalézt, vrátíme chybový výsledek
        return {-1, -1};
    }

    void initializeSubmapData(const std::vector<std::vector<int>>& submap,
                              std::unordered_set<int>& submap_set,
                              std::unordered_map<int, std::pair<int, int>>& global_to_local) {
        std::cout << "Submap content (global positions):" << std::endl;
        for (size_t x = 0; x < submap.size(); ++x) {
            for (size_t y = 0; y < submap[x].size(); ++y) {
                int global_pos = submap[x][y];
                std::cout << global_pos << " ";
                if (global_pos != -1) {
                    submap_set.insert(global_pos);
                    global_to_local[global_pos] = { static_cast<int>(x), static_cast<int>(y) };
                }
            }
            std::cout << std::endl;
        }
    }

    std::vector<std::vector<int>> generateMapRepresentation(
            const std::vector<std::vector<int>>& submap,
            const std::vector<int>& agents_in_submap,
            int problematic_timestep,
            const Instance& instance,
            const std::vector<Agent>& agents) {

        std::cout << "Map content with agents and obstacles:" << std::endl;
        std::vector<std::vector<int>> map(submap.size(), std::vector<int>(submap[0].size(), 1));
        for (size_t x = 0; x < submap.size(); ++x) {
            for (size_t y = 0; y < submap[0].size(); ++y) {
                int global_pos = submap[x][y];
                if (global_pos == -1)
                    std::cout << "X ";
                else if (instance.isObstacle(global_pos)) {
                    std::cout << "X ";
                    map[x][y] = -1;
                }
                else {
                    bool is_agent = false;
                    for (int agent : agents_in_submap) {
                        // Pozn.: agenti mohou být v submapě v různých časech, zde je použita hodnota problematic_timestep
                        if (agents[agent].path[problematic_timestep].location == global_pos) {
                            std::cout << "A ";
                            is_agent = true;
                            break;
                        }
                    }
                    if (!is_agent)
                        std::cout << ". ";
                }
            }
            std::cout << std::endl;
        }

        cout << "[DEBUG] obsah struktury map, kterou předáváme SAT solveru: " << endl;

        for (auto a : map){
            for (auto b : a){
                cout << b << "  ";
            }
            cout <<endl;
        }

        return map;
    }

    std::vector<int> getAgentsToReplan(
            const std::vector<int>& agents_in_submap,
            const std::unordered_set<int>& submap_set,
            int T_sync,
            const std::vector<Agent>& agents) {

        std::vector<int> agents_to_replan;
        std::cout << "\n[INFO] Identifikace agentů v submapě pro přeplánování "
                  << "(kdy key_agent je tam v čase " << T_sync << "):\n";

        for (int agent : agents_in_submap) {
            const auto& path = agents[agent].path;
            if (path.empty()) continue;
            int n = (int) path.size();
            std::vector<std::pair<int,int>> intervals;
            int i = 0;
            // Najdeme všechny souvislé intervaly, kdy je agent v submapě.
            while (i < n) {
                // Přeskočíme časy, kdy agent není v submapě.
                while (i < n && submap_set.find(path[i].location) == submap_set.end()) {
                    i++;
                }
                if (i >= n)
                    break;
                int t_min = i;
                // Dokud je agent v submapě, posunujeme i.
                while (i < n && submap_set.find(path[i].location) != submap_set.end()) {
                    i++;
                }
                int t_max = i - 1;
                intervals.push_back({t_min, t_max});
            }

            // Debug: vypiš všechny nalezené intervaly pro daného agenta.
            std::cout << "[DEBUG] Agent " << agent << " nalezené intervaly v submapě:";
            for (const auto& interval : intervals) {
                std::cout << " [" << interval.first << ".." << interval.second << "]";
            }
            std::cout << std::endl;

            // Vybereme interval, který obsahuje T_sync.
            bool found_interval = false;
            for (const auto& interval : intervals) {
                if (interval.first <= T_sync && T_sync <= interval.second) {
                    agents_to_replan.push_back(agent);
                    std::cout << "  Agent " << agent << " (interval v submapě: ["
                              << interval.first << ".." << interval.second << "]) obsahuje čas " << T_sync
                              << ", přidán k přeplánování.\n";
                    found_interval = true;
                    break;
                }
            }
            if (!found_interval) {
                std::cout << "[DEBUG] Agent " << agent << " není v submapě v čase " << T_sync
                          << " (intervaly: ";
                for (const auto& interval : intervals) {
                    std::cout << "[" << interval.first << ".." << interval.second << "] ";
                }
                std::cout << ")\n";
            }
        }
        if (agents_to_replan.empty())
            std::cout << "[INFO] Žádný agent nebyl v submapě ve stejnou chvíli (čas " << T_sync << ").\n";
        return agents_to_replan;
    }

    std::unordered_map<int, std::vector<std::pair<int,int>>> findLocalPaths(
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            const std::unordered_set<int>& submap_set,
            const std::unordered_map<int, std::pair<int,int>>& global_to_local,
            int T_sync,
            const std::vector<Agent>& agents) {

        // pro každého agent vytvoříme sekvenci (sx, sy) lokálních souřadnic
        // pokud agent není v submapě v čase T_sync, vynecháme ho
        // jakmile agent submapu opustí, končíme
        std::unordered_map<int, std::vector<std::pair<int,int>>> local_paths;
        std::cout << "\n[INFO] Tvorba lokálních cest (sx, sy) v submapě pro T_sync = " << T_sync << std::endl;
        for (int agent : agents_to_replan) {
            // 1) Zajistíme, že agent má definovanou pozici v T_sync
            if ((size_t)T_sync >= agents[agent].path.size()) {
                std::cout << "[WARN] Agent " << agent << " nemá definovanou pozici v čase T_sync=" << T_sync << ". Přeskakuji.\n";
                continue;
            }
            // 2) Ověříme, že je agent v submapě v T_sync
            int loc_at_Tsync = agents[agent].path[T_sync].location;
            if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
                std::cout << "[WARN] Agent " << agent << " není v submapě v čase " << T_sync << ". Přeskakuji.\n";
                continue;
            }
            // 3) Najdeme poslední čas, dokdy agent v submapě zůstává
            int last_time_in_submap = -1;
            for (int t = T_sync; t < (int)agents[agent].path.size(); t++) {
                int glob_loc = agents[agent].path[t].location;
                if (submap_set.find(glob_loc) != submap_set.end())
                    last_time_in_submap = t;
                else break;
            }
            if (last_time_in_submap == -1) {
                // Teoreticky by se to nemělo stát,
                // protože loc_at_Tsync je v submapě
                std::cout << "[WARN] Agent " << agent << " v submapě vlastně není? (podivné)\n";
                continue;
            }
            // 4) Postavíme reálnou lokální cestu v (sx, sy)
            //    od T_sync do last_time_in_submap
            std::vector<std::pair<int,int>> path_local;
            for (int t = T_sync; t <= last_time_in_submap; t++) {
                int glob_loc = agents[agent].path[t].location;
                auto it = global_to_local.find(glob_loc);
                if (it == global_to_local.end()) {
                    std::cout << "[ERROR] Agent " << agent << " v case " << t
                              << " je mimo submapu, ale last_time_in_submap=" << last_time_in_submap << std::endl;
                    break;
                }
                // (sx, sy) = lokální souřadnice v submapě
                int sx = it->second.first;
                int sy = it->second.second;
                path_local.emplace_back(sx, sy);
            }
            // 5) uložíme do mapy
            local_paths[agent] = path_local;
            std::cout << "  Agent " << agent << " (globální cesty od T=" << T_sync << " do " << last_time_in_submap << ") má lokální dráhu: ";
            for (auto& p : path_local)
                std::cout << "(" << p.first << "," << p.second << ") ";
            std::cout << std::endl;
        }
        return local_paths;
    }

    // funkce pro vypsání cesty (jednoduchý textový výpis)
    void printPathDetails(const std::vector<PathEntry>& path, int T_sync, int old_local_length, const std::vector<std::vector<int>>& submap, const std::vector<std::vector<int>>& map) {
        std::cout << "[DEBUG] Vykreslení PŮVODNÍ cesty:" << std::endl;
        // Vypiš prefix (cesta před vstupem do submapy)
        std::cout << "  Prefix (před submapou): ";
        for (int i = 0; i < T_sync && i < (int)path.size(); i++) {
            std::cout << path[i].location << " ";
        }
        std::cout << std::endl;

        // Vypiš lokální část (v submapě)
        std::cout << "  Lokální cesta (v submapě): ";
        for (int i = T_sync; i < T_sync + old_local_length && i < (int)path.size(); i++) {
            // Dekódujeme lokální pozici podle mapy
            //auto coords = decodeLocalID(path[i].location, map);
            std::cout << path[i].location << " ";
        }
        std::cout << std::endl;

        // Vypiš suffix (po opuštění submapy)
        std::cout << "  Suffix (po submapě): ";
        for (size_t i = T_sync + old_local_length; i < path.size(); i++) {
            std::cout << path[i].location << " ";
        }
        std::cout << std::endl;
    }

    bool solveWithSAT(
            std::vector<std::vector<int>>& map,
            const std::unordered_map<int, std::vector<std::pair<int,int>>>& local_paths,
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            int T_sync,
            std::vector<Agent>& agents) {

        std::cout << "\n[DEBUG] Kontrola vstupních dat pro SAT solver:\n";
        std::cout << "  - Počet agentů k přeplánování: " << agents_to_replan.size() << std::endl;

        std::vector<std::pair<int,int>> start_positions;
        std::vector<std::pair<int,int>> goal_positions;
        std::map<int,int> original_local_lengths;

        for (int agent : agents_to_replan) {
            auto it = local_paths.find(agent);
            if (it == local_paths.end() || it->second.empty()) {
                std::cout << "[WARN] Agent " << agent << " nemá local_path => vynecháváme.\n";
                continue;
            }
            const auto& path = it->second;
            start_positions.push_back(path.front());
            goal_positions.push_back(path.back());
            original_local_lengths[agent] = (int) path.size();

            std::cout << "[DEBUG] Agent " << agent << " má původní lokální dráhu délky: "
                      << path.size() << " => Start (" << path.front().first << "," << path.front().second
                      << "), Goal (" << path.back().first << "," << path.back().second << ")\n";
        }

        cout << "Startovní pozice: ";
        for (auto s : start_positions)
            cout << "(" << s.first << ", " << s.second << "), ";
        cout << endl;
        cout << "CÍLOVÉ pozice: ";
        for (auto s : goal_positions)
            cout << "(" << s.first << ", " << s.second << "), ";
        cout << endl;

        auto* inst = new _MAPFSAT_Instance(map, start_positions, goal_positions);
        auto* solver = new _MAPFSAT_DisappearAtGoal();
        auto* log = new _MAPFSAT_Logger(inst, "disappear_at_goal", 2);

        std::cout << "SAT instance and solver created.\n";

        solver->SetData(inst, log, 300, "", false, true);
        inst->SetAgents((int)start_positions.size());
        log->NewInstance((int)start_positions.size());

        int result = solver->Solve((int)start_positions.size(), 0, true, true);
        std::cout << "Solver returned: " << result << std::endl;

        if (result != 0) {
            std::cout << "SAT solver failed.\n";
            return false;
        }

        vector<vector<int>> plan = solver->GetPlan();

        for (size_t a = 0; a < plan.size(); ++a)
        {
            auto& local_path = plan[a];
            if (local_path.empty()) continue;

            /* Souřadnice cíle, který jsme solveru zadali
               (uložené při konstrukci start_positions / goal_positions). */
            const auto goal_coord = goal_positions[a];        // pair<int,int>

            // 1) Najdi PRVNÍ výskyt cíle v cestě
            size_t cut_pos = local_path.size();               // výchozí: nic nestříháme
            for (size_t t = 0; t < local_path.size(); ++t)
            {
                if (decodeLocalID(local_path[t], map) == goal_coord)
                {
                    cut_pos = t + 1;                          // necháme vrchol cíle
                    break;
                }
            }

            // 2) Odřízni vše za prvním dosažením cíle
            if (cut_pos < local_path.size())
                local_path.erase(local_path.begin() + cut_pos, local_path.end());

            // 3) Smaž případná opakování cíle na úplném konci (… 6 6 6)
            while (local_path.size() > 1 &&
                   local_path.back() == local_path[local_path.size() - 2])
            {
                local_path.pop_back();
            }
        }

        /* OD KONCE
        int i = 0;
        for (auto& path_for_agent : plan) {
            if (path_for_agent.empty()) continue;
            // Smaž opakované indexy na úplném konci (dokud se opakují v cíli):
            while (path_for_agent.size() > 1 && path_for_agent.back() == path_for_agent[path_for_agent.size() - 2]) {
                path_for_agent.pop_back();
            }
        }*/

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
            vector<PathEntry> new_local_path; // pro uložení nové lokální cesty zvlášť
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
                new_local_path.push_back(PathEntry(global_id));
                updated_path.push_back(PathEntry(global_id));
            }

            // (3) Připojíme suffix. Původní suffix začínal na indexu T_sync + old_local_length
            // (tj. tam končila stará lokální část). V nové cestě jsme skončili na indexu T_sync + new_local_length - 1
            // => suffix nalepíme od původního suffix_start dále.
            int old_suffix_start = T_sync + old_local_length;
            if (old_suffix_start < (int)agents[agent_id].path.size())
            {
                for (int t = old_suffix_start; t < (int)agents[agent_id].path.size(); t++)
                    updated_path.push_back(PathEntry(agents[agent_id].path[t].location));
            }

            // (4) Volitelná kontrola navázání prefix->lokální start
            if (T_sync > 0 && (size_t)T_sync < updated_path.size()) {
                int prefix_last = agents[agent_id].path[T_sync].location;
                int local_first = updated_path[T_sync].location;
                if (prefix_last != local_first) {
                    cout << "[WARN] agent " << agent_id
                         << " prefix->lokální navaznost se liší: "
                         << prefix_last << " != " << local_first << endl;
                }
            }

            if (T_sync + old_local_length < agents[agent_id].path.size() &&
                T_sync + new_local_length <= updated_path.size()) {
                // Lokální cíl nové lokální části je na indexu T_sync + new_local_length - 1
                int local_goal = updated_path[T_sync + new_local_length].location;
                // První prvek původního sufixu je na indexu T_sync + old_local_length (tj. kde stará lokální část končila)
                int suffix_first = agents[agent_id].path[T_sync + old_local_length].location;
                if (local_goal != suffix_first) {
                    cout << "[WARN] agent " << agent_id
                         << " lokální cíl->suffix návaznost se liší: "
                         << local_goal << " != " << suffix_first << endl;
                }
            }

            // zkontrolujeme, zda nová lokální cesta vrací původní cíl
            if (!new_local_path.empty()) {
                pair<int, int> new_local_goal = decodeLocalID(plan[a].back(), map);
                pair<int, int> original_goal = goal_positions[a];
                if (new_local_goal != original_goal) {
                    cout << "[ERROR] agent " << agent_id
                         << " nová lokální cesta končí na ("
                         << new_local_goal.first << "," << new_local_goal.second
                         << ") místo původního cíle (" << original_goal.first << ","
                         << original_goal.second << ").\n";
                    return false;
                }
            }
            else cout << "[WARN] agent " << agent_id << " nemá novou lokální cestu.\n";

            // Vykreslíme celou cestu agenta před přeplánováním, lokální a suffix:
            cout << "[DEBUG] Kompletní cesta agenta " << agent_id << ":" << endl;
            cout << "  Původní: ";
            for (auto& entry : agents[agent_id].path)
                cout << entry.location << " ";
            cout << endl;
            cout << "  Nová:     ";
            for (auto& entry : updated_path)
                cout << entry.location << " ";
            cout << endl;

            // Volitelně – voláme funkci, která vykreslí detailní informace:
            printPathDetails(updated_path, T_sync, old_local_length, submap, map);

            cout << "[INFO] Původní délka cesty agenta " << agent_id
                 << " je: " <<  agents[agent_id].path.size() << endl;

            // (5) Uložíme hotovou cestu
            agents[agent_id].path = updated_path;
            cout << "(SATUtils.cpp) nová cesta v agents[a].path agenta " << a << ": ";
            for (auto loc : agents[agent_id].path)
                cout << loc.location << ", ";
            cout << endl;

            cout << "[INFO] Cesta pro agenta " << agent_id
                 << " aktualizována, výsledná délka: "
                 << updated_path.size() << endl;
        }

        std::cout << "Paths successfully updated.\n";
        delete solver;
        delete log;
        delete inst;
        return true;
    }
}