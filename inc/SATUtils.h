//
// Created by Jaroslav Kusák on 23.03.2025.
//

#ifndef MAPF_LNS2_SATUTILS_H
#define MAPF_LNS2_SATUTILS_H


#ifndef SATUTILS_H
#define SATUTILS_H

#include <vector>
#include <utility>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <iostream>

// Předdefinice tříd – upravte dle vaší architektury
class Instance;
class PathTable;
class Agent;
class PathEntry;
class _MAPFSAT_Instance;
class _MAPFSAT_DisappearAtGoal;
class _MAPFSAT_Logger;

#include "../include/MAPF.hpp"

// V našem návrhu předpokládáme, že pro SAT solver budeme potřebovat
// také přístup k instanci a kolekci agentů. Pokud jsou tyto objekty u vás
// členy třídy LNS či InitLNS, bude třeba je předávat jako další parametry.

namespace SATUtils {

    // Vrací (sx, sy) odpovídající local_id v pořadí volných buněk v dané mapě
    std::pair<int, int> decodeLocalID(int local_id, const std::vector<std::vector<int>>& map);

    // Vytvoří submapu kolem daného agenta (na základě jeho globalID pozice) a zároveň vrátí vektor agentů,
    // kteří se v dané oblasti nacházejí. Parametry map_width a map_height definují rozměry globální mapy.
    std::pair<std::vector<std::vector<int>>, std::vector<int>> getSubmapAndAgents(
            int agent_id, int submap_size, int agent_location,
            int map_width, int map_height, const PathTable& path_table);

    // Projde zadanou submapu a naplní množinu submap_set (globalIDs platných buněk) a mapování global_to_local
    void initializeSubmapData(const std::vector<std::vector<int>>& submap,
                              std::unordered_set<int>& submap_set,
                              std::unordered_map<int, std::pair<int, int>>& global_to_local);

    // Vytvoří 2D reprezentaci mapy s ohledem na překážky a umístění agentů.
    // Funkce využívá metodu instance.isObstacle() a informace z kolekce agentů.
    std::vector<std::vector<int>> generateMapRepresentation(
            const std::vector<std::vector<int>>& submap,
            const std::vector<int>& agents_in_submap,
            int problematic_timestep,
            const Instance& instance,
            const std::vector<Agent>& agents);

    // Vybere agenty, kteří mají v submapě svůj první souvislý interval obsahující zadaný problematic_timestep.
    std::vector<int> getAgentsToReplan(
            const std::vector<int>& agents_in_submap,
            const std::unordered_set<int>& submap_set,
            int problematic_timestep,
            const std::vector<Agent>& agents);

    // Vytvoří lokální cesty pro agenty – vrací mapu, kde pro každého agenta (z agents_to_replan)
    // je vektor dvojic (sx, sy) odpovídající jeho pozicím v submapě od T_sync do okamžiku, kdy opustí submapu.
    std::unordered_map<int, std::vector<std::pair<int,int>>> findLocalPaths(
            const std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            const std::unordered_set<int>& submap_set,
            const std::unordered_map<int, std::pair<int,int>>& global_to_local,
            int T_sync,
            const std::vector<Agent>& agents);

    // Zavolá SAT solver – vytvoří SAT instanci, spustí solver, dekóduje výsledky a aktualizuje cesty.
    // Vrací true, pokud se podařilo najít nové řešení.
    bool solveWithSAT(
            std::vector<std::vector<int>>& map,
            const std::unordered_map<int, std::vector<std::pair<int,int>>>& local_paths,
            std::vector<int>& agents_to_replan,
            const std::vector<std::vector<int>>& submap,
            int T_sync,
            std::vector<Agent>& agents);
}

#endif // SATUTILS_H


#endif //MAPF_LNS2_SATUTILS_H
