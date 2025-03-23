#include "SATUtils.h"

// V případě potřeby přidejte další includy – např. pro SAT solver, Instance, Agent, PathTable apod.
#include <iostream>
#include <algorithm>

// Implementace funkcí v rámci namespace SATUtils
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
                    std::cout << ". ";
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
        return map;
    }

    std::vector<int> getAgentsToReplan(
            const std::vector<int>& agents_in_submap,
            const std::unordered_set<int>& submap_set,
            int problematic_timestep,
            const std::vector<Agent>& agents) {

        std::vector<int> agents_to_replan;
        std::cout << "\n[INFO] Identifikace agentů v submapě pro přeplánování "
                  << "(kdy key_agent je tam v čase " << problematic_timestep << "):\n";

        for (int agent : agents_in_submap) {
            const auto& path = agents[agent].path;
            if (path.empty()) continue;
            int t_min = -1, t_max = -1;
            for (int t = 0; t < (int)path.size(); t++) {
                int loc = path[t].location;
                bool inSubmap = (submap_set.find(loc) != submap_set.end());
                if (t_min < 0) {
                    if (inSubmap) {
                        t_min = t;
                        t_max = t;
                    }
                }
                else {
                    if (inSubmap)
                        t_max = t;
                    else break;
                }
            }
            if (t_min == -1)
                continue;
            if (t_min <= problematic_timestep && problematic_timestep <= t_max) {
                agents_to_replan.push_back(agent);
                std::cout << "  Agent " << agent << " (první interval v submapě je ["
                          << t_min << ".." << t_max << "]), pokrývá i čas " << problematic_timestep
                          << ", přidán k přeplánování.\n";
            }
        }
        if (agents_to_replan.empty())
            std::cout << "[INFO] Žádný agent nebyl v submapě ve stejnou chvíli (čas "
                      << problematic_timestep << ").\n";
        return agents_to_replan;
    }

    std::unordered_map<int, std::vector<std::pair<int,int>>> findLocalPaths(
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            const std::unordered_set<int>& submap_set,
            const std::unordered_map<int, std::pair<int,int>>& global_to_local,
            int T_sync,
            const std::vector<Agent>& agents) {

        std::unordered_map<int, std::vector<std::pair<int,int>>> local_paths;
        std::cout << "\n[INFO] Tvorba lokálních cest (sx, sy) v submapě pro T_sync = " << T_sync << std::endl;
        for (int agent : agents_to_replan) {
            if ((size_t)T_sync >= agents[agent].path.size()) {
                std::cout << "[WARN] Agent " << agent << " nemá definovanou pozici v čase T_sync=" << T_sync << ". Přeskakuji.\n";
                continue;
            }
            int loc_at_Tsync = agents[agent].path[T_sync].location;
            if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
                std::cout << "[WARN] Agent " << agent << " není v submapě v čase " << T_sync << ". Přeskakuji.\n";
                continue;
            }
            int last_time_in_submap = -1;
            for (int t = T_sync; t < (int)agents[agent].path.size(); t++) {
                int glob_loc = agents[agent].path[t].location;
                if (submap_set.find(glob_loc) != submap_set.end())
                    last_time_in_submap = t;
                else break;
            }
            if (last_time_in_submap == -1) {
                std::cout << "[WARN] Agent " << agent << " v submapě vlastně není? (podivné)\n";
                continue;
            }
            std::vector<std::pair<int,int>> path_local;
            for (int t = T_sync; t <= last_time_in_submap; t++) {
                int glob_loc = agents[agent].path[t].location;
                auto it = global_to_local.find(glob_loc);
                if (it == global_to_local.end()) {
                    std::cout << "[ERROR] Agent " << agent << " v case " << t
                              << " je mimo submapu, ale last_time_in_submap=" << last_time_in_submap << std::endl;
                    break;
                }
                int sx = it->second.first;
                int sy = it->second.second;
                path_local.emplace_back(sx, sy);
            }
            local_paths[agent] = path_local;
            std::cout << "  Agent " << agent << " (globální cesty od T=" << T_sync << " do " << last_time_in_submap << ") má lokální dráhu: ";
            for (auto& p : path_local)
                std::cout << "(" << p.first << "," << p.second << ") ";
            std::cout << std::endl;
        }
        return local_paths;
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

        auto inst = std::make_unique<_MAPFSAT_Instance>(map, start_positions, goal_positions);
        auto solver = std::make_unique<_MAPFSAT_DisappearAtGoal>();
        auto log = std::make_unique<_MAPFSAT_Logger>(inst.get(), "disappear_at_goal", 2);

        std::cout << "SAT instance and solver created.\n";

        solver->SetData(inst.get(), log.get(), 300, "", false, true);
        inst->SetAgents((int)start_positions.size());
        log->NewInstance((int)start_positions.size());

        int result = solver->Solve((int)start_positions.size(), 0, true);
        std::cout << "Solver returned: " << result << std::endl;

        if (result != 0) {
            std::cout << "SAT solver failed.\n";
            return false;
        }

        std::vector<std::vector<int>> plan = solver->GetPlan();

        for (auto& path_for_agent : plan) {
            if (path_for_agent.empty()) continue;
            while (path_for_agent.size() > 1 &&
                   path_for_agent.back() == path_for_agent[path_for_agent.size() - 2]) {
                path_for_agent.pop_back();
            }
        }

        for (size_t a = 0; a < plan.size(); ++a) {
            std::cout << "[DEBUG] Agent (index) " << agents_to_replan[a]
                      << " | Nová lokální cesta (submap idx): ";
            for (auto lid : plan[a])
                std::cout << lid << " ";
            std::cout << std::endl;
        }

        for (size_t a = 0; a < plan.size(); ++a) {
            int agent_id = agents_to_replan[a];
            int old_local_length = original_local_lengths[agent_id];
            int new_local_length = (int) plan[a].size();

            std::cout << "[INFO] Aktualizace cesty pro agenta " << agent_id
                      << " | Původní lokální délka: " << old_local_length
                      << " | Nová lokální délka: " << new_local_length << std::endl;

            std::vector<PathEntry> updated_path(
                    agents[agent_id].path.begin(),
                    agents[agent_id].path.begin() + T_sync
            );

            for (int t = 0; t < new_local_length; t++) {
                int local_id = plan[a][t];
                std::pair<int,int> coords = decodeLocalID(local_id, map);
                int sx = coords.first, sy = coords.second;
                if (sx == -1 || sy == -1) {
                    std::cout << "[ERROR] agent " << agent_id
                              << " local_id=" << local_id << " nelze dekódovat do platných souřadnic.\n";
                    continue;
                }
                int global_id = submap[sx][sy];
                std::cout << "[DEBUG] agent " << agent_id << " t=" << t
                          << " => decoded (sx,sy)=(" << sx << "," << sy
                          << ") => global_id=" << global_id << std::endl;
                updated_path.push_back(PathEntry(global_id));
            }

            int old_suffix_start = T_sync + old_local_length;
            if (old_suffix_start < (int)agents[agent_id].path.size()) {
                for (int t = old_suffix_start; t < (int)agents[agent_id].path.size(); t++) {
                    updated_path.push_back(PathEntry(agents[agent_id].path[t].location));
                }
            }

            if (T_sync > 0 && (size_t)T_sync < updated_path.size()) {
                int prefix_last = agents[agent_id].path[T_sync].location;
                int local_first = updated_path[T_sync].location;
                if (prefix_last != local_first) {
                    std::cout << "[WARN] agent " << agent_id
                              << " prefix->lokální navaznost se liší: "
                              << prefix_last << " != " << local_first << std::endl;
                }
            }

            std::cout << "[INFO] Původní délka cesty agenta " << agent_id
                      << " je: " << agents[agent_id].path.size() << std::endl;
            agents[agent_id].path = updated_path;
            std::cout << "[INFO] Cesta pro agenta " << agent_id
                      << " aktualizována, výsledná délka: " << updated_path.size() << std::endl;
        }

        std::cout << "Paths successfully updated.\n";
        return true;
    }
}