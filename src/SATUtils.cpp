#include "SATUtils.h"
#include "Log.h"

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
        SAT_DBG("Submap content (global positions):");
        for (size_t x = 0; x < submap.size(); ++x) {
            for (size_t y = 0; y < submap[x].size(); ++y) {
                int global_pos = submap[x][y];
                if (global_pos != -1) {
                    submap_set.insert(global_pos);
                    global_to_local[global_pos] = { static_cast<int>(x), static_cast<int>(y) };
                }
            }
        }
        if (satlog::debug_enabled) {
            for (size_t x = 0; x < submap.size(); ++x) {
                std::stringstream ss;
                for (size_t y = 0; y < submap[x].size(); ++y) {
                    int global_pos = submap[x][y];
                    ss << global_pos << " ";
                }
                SAT_DBG(ss.str());
            }
        }
    }

    // Generates a map representation with agents and obstacles for the SAT solver.
    std::vector<std::vector<int>> generateMapRepresentation(
            const std::vector<std::vector<int>>& submap,
            const std::vector<int>& agents_in_submap,
            int problematic_timestep,
            const Instance& instance,
            const std::vector<Agent>& agents) {

        SAT_DBG("Map content with agents and obstacles:");
        std::vector<std::vector<int>> map(submap.size(), std::vector<int>(submap[0].size(), 1));
        if (satlog::debug_enabled) {
            for (size_t x = 0; x < submap.size(); ++x) {
                std::stringstream ss;
                for (size_t y = 0; y < submap[0].size(); ++y) {
                    int global_pos = submap[x][y];
                    if (global_pos == -1)
                        ss << "X ";
                    else if (instance.isObstacle(global_pos)) {
                        ss << "X ";
                    }
                    else {
                        bool is_agent = false;
                        for (int agent : agents_in_submap) {
                            // Note: agents can occupy the sub‑map at different timesteps; we reference problematic_timestep here
                            if (agents[agent].path[problematic_timestep].location == global_pos) {
                                ss << "A ";
                                is_agent = true;
                                break;
                            }
                        }
                        if (!is_agent)
                            ss << ". ";
                    }
                }
                SAT_DBG(ss.str());
            }
        }
        for (size_t x = 0; x < submap.size(); ++x) {
            for (size_t y = 0; y < submap[0].size(); ++y) {
                int global_pos = submap[x][y];
                if (global_pos == -1) {
                    map[x][y] = 1; // keep as 1 (or could be -1, but follows previous logic)
                }
                else if (instance.isObstacle(global_pos)) {
                    map[x][y] = -1;
                }
                // else: leave as 1 (free)
            }
        }

        SAT_DBG("content of the map structure that we pass to the SAT solver:");
        if (satlog::debug_enabled) {
            for (auto& row : map) {
                std::stringstream ss;
                for (auto cell : row) ss << cell << "  ";
                SAT_DBG(ss.str());
            }
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
        SAT_DBG("Identifying agents in a submap for rescheduling (when key_agent is there at the time " << T_sync << "):");

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
            SAT_DBG("Agent " << agent << " found intervals in the sub-map:");
            if (satlog::debug_enabled) {
                std::stringstream ss;
                for (const auto& interval : intervals)
                    ss << " [" << interval.first << ".." << interval.second << "]";
                SAT_DBG(ss.str());
            }

            // Select the interval that contains T_sync.
            bool found_interval = false;
            for (const auto& interval : intervals) {
                if (interval.first <= T_sync && T_sync <= interval.second) {
                    agents_to_replan.push_back(agent);
                    SAT_DBG("Agent " << agent << " (interval in submap: ["
                              << interval.first << ".." << interval.second << "]) contains time " << T_sync
                              << ", added to replan.");
                    found_interval = true;
                    break;
                }
            }
            if (!found_interval) {
                SAT_DBG("Agent " << agent << " is not in sub-map at time " << T_sync << " (intervals: ");
                if (satlog::debug_enabled) {
                    std::stringstream ss;
                    for (const auto& interval : intervals) {
                        ss << "[" << interval.first << ".." << interval.second << "] ";
                    }
                    ss << ")";
                    SAT_DBG(ss.str());
                }
            }
        }
        if (agents_to_replan.empty())
            SAT_DBG("No agent was in the submap at the same time (time " << T_sync << ").");
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
        SAT_DBG("Creating local paths (sx, sy) in the submap for T_sync = " << T_sync);
        for (int agent : agents_to_replan) {
            // Ensure the agent has a defined position at T_sync
            if ((size_t)T_sync >= agents[agent].path.size()) {
                SAT_DBG("Agent " << agent << " has no defined position in time T_sync=" << T_sync << ". Skipping.");
                continue;
            }
            // Verify the agent is in the sub‑map at T_sync
            int loc_at_Tsync = agents[agent].path[T_sync].location;
            if (submap_set.find(loc_at_Tsync) == submap_set.end()) {
                SAT_DBG("Agent " << agent << " not in submap in time " << T_sync << ". Skipping.");
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
                SAT_DBG("Agent " << agent << " It's actually not in the submap? (strange)");
                continue;
            }
            // Build the true local path in (sx, sy)
            // from T_sync to last_time_in_submap
            std::vector<std::pair<int,int>> path_local;
            for (int t = T_sync; t <= last_time_in_submap; t++) {
                int glob_loc = agents[agent].path[t].location;
                auto it = global_to_local.find(glob_loc);
                if (it == global_to_local.end()) {
                    SAT_DBG("Agent " << agent << " at time " << t
                              << " is out of submap but last_time_in_submap=" << last_time_in_submap);
                    break;
                }
                // (sx, sy) = local coordinates in the sub‑map
                int sx = it->second.first;
                int sy = it->second.second;
                path_local.emplace_back(sx, sy);
            }
            // store into the resulting map
            local_paths[agent] = path_local;
            {
                std::stringstream ss;
                ss << "  Agent " << agent << " (global paths from T=" << T_sync << " to " << last_time_in_submap << ") has local path: ";
                for (auto& p : path_local)
                    ss << "(" << p.first << "," << p.second << ") ";
                SAT_DBG(ss.str());
            }
        }
        return local_paths;
    }

    // Prints detailed information about an agent's path (prefix, local, and suffix segments).
    void printPathDetails(const std::vector<PathEntry>& path, int T_sync, int old_local_length, const std::vector<std::vector<int>>& submap, const std::vector<std::vector<int>>& map) {
        SAT_DBG("Rendering the ORIGINAL path:");
        // Print prefix (path before entering the sub‑map)
        {
            std::stringstream ss_prefix;
            ss_prefix << "  Prefix (before submap): ";
            for (int i = 0; i < T_sync && i < (int)path.size(); i++) {
                ss_prefix << path[i].location << " ";
            }
            SAT_DBG(ss_prefix.str());
        }
        {
            std::stringstream ss_local;
            ss_local << "  Local path (in submap): ";
            for (int i = T_sync; i < T_sync + old_local_length && i < (int)path.size(); i++) {
                ss_local << path[i].location << " ";
            }
            SAT_DBG(ss_local.str());
        }
        {
            std::stringstream ss_suffix;
            ss_suffix << "  Suffix (out of submap): ";
            for (size_t i = T_sync + old_local_length; i < path.size(); i++) {
                ss_suffix << path[i].location << " ";
            }
            SAT_DBG(ss_suffix.str());
        }
    }

    // Solves the local MAPF problem in the sub-map using a SAT-based solver and updates agent paths.
    bool solveWithSAT(
            std::vector<std::vector<int>>& map,
            const std::unordered_map<int, std::vector<std::pair<int,int>>>& local_paths,
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            int T_sync,
            std::vector<Agent>& agents) {

        SAT_DBG("Checking input data for SAT solver:");
        SAT_DBG("  - Number of agents to replan: " << agents_to_replan.size());

        std::vector<std::pair<int,int>> start_positions;
        std::vector<std::pair<int,int>> goal_positions;
        std::map<int,int> original_local_lengths;

        for (int agent : agents_to_replan) {
            auto it = local_paths.find(agent);
            if (it == local_paths.end() || it->second.empty()) {
                SAT_DBG("Agent " << agent << " doesn't have local_path => skipping.");
                continue;
            }
            const auto& path = it->second;
            start_positions.push_back(path.front());
            goal_positions.push_back(path.back());
            original_local_lengths[agent] = (int) path.size();

            SAT_DBG("Agent " << agent << " has the original local path length: "
                      << path.size() << " => Start (" << path.front().first << "," << path.front().second
                      << "), Goal (" << path.back().first << "," << path.back().second << ")");
        }

        {
            std::stringstream ss;
            ss << "Start position: ";
            for (auto s : start_positions)
                ss << "(" << s.first << ", " << s.second << "), ";
            SAT_DBG(ss.str());
        }
        {
            std::stringstream ss;
            ss << "Goal position: ";
            for (auto s : goal_positions)
                ss << "(" << s.first << ", " << s.second << "), ";
            SAT_DBG(ss.str());
        }

        auto inst = std::make_unique<_MAPFSAT_Instance>(map, start_positions, goal_positions);
        auto solver = std::make_unique<_MAPFSAT_DisappearAtGoal>();
        auto log = std::make_unique<_MAPFSAT_Logger>(inst.get(), "disappear_at_goal", 2);

        SAT_DBG("SAT instance and solver created.");

        solver->SetData(inst.get(), log.get(), 300, "", false, true);
        inst->SetAgents((int)start_positions.size());
        log->NewInstance((int)start_positions.size());

        int result = solver->Solve((int)start_positions.size(), 0, true, true);
        SAT_STAT("Solver returned: " << result);

        if (result != 0) {
            SAT_DBG("SAT solver failed.");
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
            SAT_DBG("Agent (index) " << agents_to_replan[a] << " | New local path (submap idx):");
            if (satlog::debug_enabled) {
                std::stringstream ss;
                for (auto lid : plan[a]) {
                    ss << lid << " ";
                }
                SAT_DBG(ss.str());
            }
        }

        // Update updated_path for each agent (prefix + new local path + suffix)
        for (size_t a = 0; a < plan.size(); ++a) {
            int agent_id = agents_to_replan[a];

            // (A) Original local segment length
            int old_local_length = original_local_lengths[agent_id];
            // (B) New local segment length
            int new_local_length = (int) plan[a].size();

            SAT_DBG("Update the agent path " << agent_id
                 << " | Original local length: " << old_local_length
                 << " | New local length: " << new_local_length);

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
                    SAT_DBG("agent " << agent_id
                         << " local_id=" << local_id
                         << " cannot be decoded into valid coordinates.");
                    continue;
                }
                // Global ID from the sub‑map:
                int global_id = submap[sx][sy];
            SAT_DBG("agent " << agent_id << " t=" << t
                 << " => decoded (sx,sy)=(" << sx << "," << sy
                 << ") => global_id=" << global_id);
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
            else SAT_DBG("agent " << agent_id << " does not have a new local path.");

            // Print the agent's entire path before replanning, the local segment, and the suffix:
            SAT_DBG("Complete path of agent " << agent_id << ":");
            {
                std::stringstream ss;
                ss << "  Original: ";
                for (auto& entry : agents[agent_id].path)
                    ss << entry.location << " ";
                SAT_DBG(ss.str());
            }
            {
                std::stringstream ss;
                ss << "  New:     ";
                for (auto& entry : updated_path)
                    ss << entry.location << " ";
                SAT_DBG(ss.str());
            }

            // optionally call the helper that prints detailed information:
            printPathDetails(updated_path, T_sync, old_local_length, submap, map);

            SAT_DBG("Original agent path length " << agent_id
                 << " je: " <<  agents[agent_id].path.size());

            // save the final path
            agents[agent_id].path = updated_path;
            {
                std::stringstream ss;
                ss << "(SATUtils.cpp) new path in agents[a].path of the agent " << a << ": ";
                for (auto loc : agents[agent_id].path)
                    ss << loc.location << ", ";
                SAT_DBG(ss.str());
            }

            SAT_DBG("Path for agent " << agent_id
                 << " updated, resulting length: "
                 << updated_path.size());
        }

        SAT_DBG("Paths successfully updated.");
        return true;
    }
}