#include "SATUtils.h"

#include <iostream>
#include <algorithm>

namespace SATUtils {
    // Decodes a local ID into (row, column) coordinates in the map for a free cell.
    std::pair<int, int> decodeLocalID(int local_id, const std::vector<std::vector<int>>& map) {
        int count = 0;
        int rows = map.size();
        int cols = map[0].size();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (map[i][j] == 1) { // if the cell is free
                    if (count == local_id)
                        return {i, j};
                    count++;
                }
            }
        }
        // If the index could not be found, return an error pair
        return {-1, -1};
    }

    // Initializes sub-map data structures for quick lookup and coordinate mapping.
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

    // Generates a map representation with agents and obstacles for the SAT solver.
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
                        // Note: agents can occupy the sub‑map at different timesteps; we reference problematic_timestep here
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

        cout << "[DEBUG] content of the map structure that we pass to the SAT solver: " << endl;

        for (auto a : map){
            for (auto b : a)
                cout << b << "  ";
            cout <<endl;
        }

        return map;
    }

    // Finds agents that need to be replanned based on their presence in the sub-map at T_sync.
    std::vector<int> getAgentsToReplan(
            const std::vector<int>& agents_in_submap,
            const std::unordered_set<int>& submap_set,
            int T_sync,
            const std::vector<Agent>& agents) {

        std::vector<int> agents_to_replan;
        std::cout << "\n[INFO] Identifying agents in a submap for rescheduling "
                  << "(when key_agent is there at the time " << T_sync << "):\n";

        for (int agent : agents_in_submap) {
            const auto& path = agents[agent].path;
            if (path.empty()) continue;
            int n = (int) path.size();
            std::vector<std::pair<int,int>> intervals;
            int i = 0;
            // Find all contiguous intervals when the agent stays inside the sub‑map.
            while (i < n) {
                // Skip timesteps when the agent is outside the sub‑map.
                while (i < n && submap_set.find(path[i].location) == submap_set.end()) {
                    i++;
                }
                if (i >= n)
                    break;
                int t_min = i;
                // Advance while the agent remains inside the sub‑map.
                while (i < n && submap_set.find(path[i].location) != submap_set.end()) {
                    i++;
                }
                int t_max = i - 1;
                intervals.push_back({t_min, t_max});
            }

            // Debug: print all intervals found for this agent:
            std::cout << "[DEBUG] Agent " << agent << " found intervals in the submap:";
            for (const auto& interval : intervals) {
                std::cout << " [" << interval.first << ".." << interval.second << "]";
            }
            std::cout << std::endl;

            // Select the interval that contains T_sync.
            bool found_interval = false;
            for (const auto& interval : intervals) {
                if (interval.first <= T_sync && T_sync <= interval.second) {
                    agents_to_replan.push_back(agent);
                    std::cout << "  Agent " << agent << " (interval in submap: ["
                              << interval.first << ".." << interval.second << "]) contains time " << T_sync
                              << ", added to replan.\n";
                    found_interval = true;
                    break;
                }
            }
            if (!found_interval) {
                std::cout << "[DEBUG] Agent " << agent << " is not in submap at time " << T_sync
                          << " (intervals: ";
                for (const auto& interval : intervals) {
                    std::cout << "[" << interval.first << ".." << interval.second << "] ";
                }
                std::cout << ")\n";
            }
        }
        if (agents_to_replan.empty())
            std::cout << "[INFO] No agent was in the submap at the same time (time " << T_sync << ").\n";
        return agents_to_replan;
    }

    // For every agent build a sequence (sx, sy) of local coordinates. Skip agents not present in the sub‑map at T_sync; stop once the agent leaves the sub‑map.
    std::unordered_map<int, std::vector<std::pair<int,int>>> findLocalPaths(
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            const std::unordered_set<int>& submap_set,
            const std::unordered_map<int, std::pair<int,int>>& global_to_local,
            int T_sync,
            const std::vector<Agent>& agents) {

        std::unordered_map<int, std::vector<std::pair<int,int>>> local_paths;
        std::cout << "\n[INFO] Creating local paths (sx, sy) in the submap for T_sync = " << T_sync << std::endl;
        for (int agent : agents_to_replan) {
            // Ensure the agent has a defined position at T_sync
            if ((size_t)T_sync >= agents[agent].path.size()) {
                std::cout << "[WARN] Agent " << agent << " has no defined position in time T_sync=" << T_sync << ". Skipping.\n";
                continue;
            }
            // Verify the agent is in the sub‑map at T_sync
            int loc_at_Tsync = agents[agent].path[T_sync].location;
            if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
                std::cout << "[WARN] Agent " << agent << " not in submap in time " << T_sync << ". Skipping.\n";
                continue;
            }
            // Find the last timestep the agent remains in the sub‑map
            int last_time_in_submap = -1;
            for (int t = T_sync; t < (int)agents[agent].path.size(); t++) {
                int glob_loc = agents[agent].path[t].location;
                if (submap_set.find(glob_loc) != submap_set.end())
                    last_time_in_submap = t;
                else break;
            }
            if (last_time_in_submap == -1) {
                // Theoretically this should not happen because loc_at_Tsync is inside the sub‑map
                std::cout << "[WARN] Agent " << agent << " It's actually not in the submap? (strange)\n";
                continue;
            }
            // Build the true local path in (sx, sy)
            // from T_sync to last_time_in_submap
            std::vector<std::pair<int,int>> path_local;
            for (int t = T_sync; t <= last_time_in_submap; t++) {
                int glob_loc = agents[agent].path[t].location;
                auto it = global_to_local.find(glob_loc);
                if (it == global_to_local.end()) {
                    std::cout << "[ERROR] Agent " << agent << " at time " << t
                              << " is out of submap but last_time_in_submap=" << last_time_in_submap << std::endl;
                    break;
                }
                // (sx, sy) = local coordinates in the sub‑map
                int sx = it->second.first;
                int sy = it->second.second;
                path_local.emplace_back(sx, sy);
            }
            // store into the resulting map
            local_paths[agent] = path_local;
            std::cout << "  Agent " << agent << " (global paths from T=" << T_sync << " to " << last_time_in_submap << ") has local path: ";
            for (auto& p : path_local)
                std::cout << "(" << p.first << "," << p.second << ") ";
            std::cout << std::endl;
        }
        return local_paths;
    }

    // Prints detailed information about an agent's path (prefix, local, and suffix segments).
    void printPathDetails(const std::vector<PathEntry>& path, int T_sync, int old_local_length, const std::vector<std::vector<int>>& submap, const std::vector<std::vector<int>>& map) {
        std::cout << "[DEBUG] Rendering the ORIGINAL path:" << std::endl;
        // Print prefix (path before entering the sub‑map)
        std::cout << "  Prefix (before submap): ";
        for (int i = 0; i < T_sync && i < (int)path.size(); i++) {
            std::cout << path[i].location << " ";
        }
        std::cout << std::endl;

        // Print local segment (inside sub‑map)
        std::cout << "  Local path (in submap): ";
        for (int i = T_sync; i < T_sync + old_local_length && i < (int)path.size(); i++) {
            // Decode the local position using the map
            //auto coords = decodeLocalID(path[i].location, map);
            std::cout << path[i].location << " ";
        }
        std::cout << std::endl;

        // Print suffix (after leaving sub‑map)
        std::cout << "  Suffix (out of submap): ";
        for (size_t i = T_sync + old_local_length; i < path.size(); i++) {
            std::cout << path[i].location << " ";
        }
        std::cout << std::endl;
    }

    // Solves the local MAPF problem in the sub-map using a SAT-based solver and updates agent paths.
    bool solveWithSAT(
            std::vector<std::vector<int>>& map,
            const std::unordered_map<int, std::vector<std::pair<int,int>>>& local_paths,
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            int T_sync,
            std::vector<Agent>& agents) {

        std::cout << "\n[DEBUG] Checking input data for SAT solver:\n";
        std::cout << "  - Number of agents to replan: " << agents_to_replan.size() << std::endl;

        std::vector<std::pair<int,int>> start_positions;
        std::vector<std::pair<int,int>> goal_positions;
        std::map<int,int> original_local_lengths;

        for (int agent : agents_to_replan) {
            auto it = local_paths.find(agent);
            if (it == local_paths.end() || it->second.empty()) {
                std::cout << "[WARN] Agent " << agent << " doesn't have local_path => skipping.\n";
                continue;
            }
            const auto& path = it->second;
            start_positions.push_back(path.front());
            goal_positions.push_back(path.back());
            original_local_lengths[agent] = (int) path.size();

            std::cout << "[DEBUG] Agent " << agent << " has the original local path length: "
                      << path.size() << " => Start (" << path.front().first << "," << path.front().second
                      << "), Goal (" << path.back().first << "," << path.back().second << ")\n";
        }

        cout << "Start position: ";
        for (auto s : start_positions)
            cout << "(" << s.first << ", " << s.second << "), ";
        cout << endl;
        cout << "Goal position: ";
        for (auto s : goal_positions)
            cout << "(" << s.first << ", " << s.second << "), ";
        cout << endl;

        auto inst = std::make_unique<_MAPFSAT_Instance>(map, start_positions, goal_positions);
        auto solver = std::make_unique<_MAPFSAT_DisappearAtGoal>();
        auto log = std::make_unique<_MAPFSAT_Logger>(inst.get(), "disappear_at_goal", 2);

        std::cout << "SAT instance and solver created.\n";

        solver->SetData(inst.get(), log.get(), 300, "", false, true);
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

            /* Coordinates of the goal supplied to the solver
               (recorded when start_positions / goal_positions were created). */
            const auto goal_coord = goal_positions[a];        // pair<int,int>

            // Find the FIRST occurrence of the goal in the path
            size_t cut_pos = local_path.size();               // default: not cutting out anything
            for (size_t t = 0; t < local_path.size(); ++t) {
                if (decodeLocalID(local_path[t], map) == goal_coord) {
                    cut_pos = t + 1;                          // we leave the final node
                    break;
                }
            }

            // Remove everything after the first arrival at the goal
            if (cut_pos < local_path.size())
                local_path.erase(local_path.begin() + cut_pos, local_path.end());

            // Trim possible repetitions of the goal at the very end (e.g., 6 6 6)
            while (local_path.size() > 1 &&
                   local_path.back() == local_path[local_path.size() - 2])
            {
                local_path.pop_back();
            }
        }

        // Debug: print new local path (in 1D submap indices).
        for (size_t a = 0; a < plan.size(); ++a) {
            cout << "[DEBUG] Agent (index) " << agents_to_replan[a]
                 << " | New local path (submap idx): ";
            for (auto lid : plan[a]) {
                cout << lid << " ";
            }
            cout << endl;
        }

        // Update updated_path for each agent (prefix + new local path + suffix)
        for (size_t a = 0; a < plan.size(); ++a) {
            int agent_id = agents_to_replan[a];

            // (A) Original local segment length
            int old_local_length = original_local_lengths[agent_id];
            // (B) New local segment length
            int new_local_length = (int) plan[a].size();

            cout << "[INFO] Update the agent path " << agent_id
                 << " | Original local length: " << old_local_length
                 << " | New local length: " << new_local_length << endl;

            // Copy the prefix up to T_sync
            vector<PathEntry> updated_path(
                    agents[agent_id].path.begin(),
                    agents[agent_id].path.begin() + T_sync
            );

            // Insert the new local segment (decoded to global IDs)
            vector<PathEntry> new_local_path; // to save a new local path separately
            for (int t = 0; t < new_local_length; t++)
            {
                int local_id = plan[a][t];
                // decodeLocalID => (sx,sy) submap
                // *** WARNING: we have to ensure correct x/y orientation and submap[x][y] vs submap[y][x]! ***
                pair<int, int> coords = decodeLocalID(local_id, map);
                int sx = coords.first;
                int sy = coords.second;
                if (sx == -1 || sy == -1) {
                    cout << "[ERROR] agent " << agent_id
                         << " local_id=" << local_id
                         << " cannot be decoded into valid coordinates.\n";
                    continue;
                }
                // Global ID from the sub‑map:
                int global_id = submap[sx][sy];
                cout << "[DEBUG] agent " << agent_id << " t=" << t
                     << " => decoded (sx,sy)=(" << sx << "," << sy
                     << ") => global_id=" << global_id << endl;
                new_local_path.push_back(PathEntry(global_id));
                updated_path.push_back(PathEntry(global_id));
            }

            // Append the suffix. The original suffix started at index T_sync + old_local_length; the new path ends at T_sync + new_local_length - 1
            int old_suffix_start = T_sync + old_local_length;
            if (old_suffix_start < (int)agents[agent_id].path.size()) {
                for (int t = old_suffix_start; t < (int)agents[agent_id].path.size(); t++)
                    updated_path.push_back(PathEntry(agents[agent_id].path[t].location));
            }

            // Optional check of prefix -> local start continuity
            if (T_sync > 0 && (size_t)T_sync < updated_path.size()) {
                int prefix_last = agents[agent_id].path[T_sync].location;
                int local_first = updated_path[T_sync].location;
                if (prefix_last != local_first) {
                    cout << "[WARN] agent " << agent_id
                         << " prefix -> local continuity varies: "
                         << prefix_last << " != " << local_first << endl;
                }
            }

            if (T_sync + old_local_length < agents[agent_id].path.size() &&
                T_sync + new_local_length <= updated_path.size()) {
                // The local goal of the new segment is at index T_sync + new_local_length - 1
                int local_goal = updated_path[T_sync + new_local_length].location;
                // The first element of the original suffix is at index T_sync + old_local_length (i.e. where the old local part ended)                int suffix_first = agents[agent_id].path[T_sync + old_local_length].location;
                int suffix_first = agents[agent_id].path[T_sync + old_local_length].location;
                if (local_goal != suffix_first) {
                    cout << "[WARN] agent " << agent_id
                         << " local goal -> suffix continuity varies: "
                         << local_goal << " != " << suffix_first << endl;
                }
            }

            // Verify that the new local path matches the original goal
            if (!new_local_path.empty()) {
                pair<int, int> new_local_goal = decodeLocalID(plan[a].back(), map);
                pair<int, int> original_goal = goal_positions[a];
                if (new_local_goal != original_goal) {
                    cout << "[ERROR] agent " << agent_id
                         << " new local path ends at ("
                         << new_local_goal.first << "," << new_local_goal.second
                         << ") instead of original goal (" << original_goal.first << ","
                         << original_goal.second << ").\n";
                    return false;
                }
            }
            else cout << "[WARN] agent " << agent_id << " does not have a new local path.\n";

            // Print the agent's entire path before replanning, the local segment, and the suffix:
            cout << "[DEBUG] Complete path of agent " << agent_id << ":" << endl;
            cout << "  Original: ";
            for (auto& entry : agents[agent_id].path)
                cout << entry.location << " ";
            cout << endl;
            cout << "  New:     ";
            for (auto& entry : updated_path)
                cout << entry.location << " ";
            cout << endl;

            // Optionally call the helper that prints detailed information:
            printPathDetails(updated_path, T_sync, old_local_length, submap, map);

            cout << "[INFO] Original agent path length " << agent_id
                 << " je: " <<  agents[agent_id].path.size() << endl;

            // Save the final path
            agents[agent_id].path = updated_path;
            cout << "(SATUtils.cpp) new path in agents[a].path of the agent " << a << ": ";
            for (auto loc : agents[agent_id].path)
                cout << loc.location << ", ";
            cout << endl;

            cout << "[INFO] Path for agent " << agent_id
                 << " updated, resulting length: "
                 << updated_path.size() << endl;
        }

        std::cout << "Paths successfully updated.\n";
        return true;
    }
}