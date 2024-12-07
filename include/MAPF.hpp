#ifndef _instance_h_INCLUDED
#define _instance_h_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <queue>
#include <algorithm>

struct _MAPFSAT_Vertex
{
	size_t x;
	size_t y;

	bool operator==(const _MAPFSAT_Vertex &rhs) const
	{
		if (rhs.x == x && rhs.y == y)
			return true;
		return false;
    }
};

struct _MAPFSAT_Agent
{
	_MAPFSAT_Vertex start;
	_MAPFSAT_Vertex goal;
};

struct _MAPFSAT_Avoid
{
	_MAPFSAT_Vertex v;
	int t;
};

class _MAPFSAT_Instance
{
public:
	/** Constructor of _MAPFSAT_Instance using files.
    *
    * Stores all relevant data about the solved instance.
    *
    * @param map_dir directory including the map.
    * @param agents_file file containing the solved scenario.
    */
    _MAPFSAT_Instance(std::string, std::string);

    /** Constructor of _MAPFSAT_Instance using data.
    *
    * Stores all relevant data about the solved instance.
    *
    * @param map_vector the map, obstacles are marked by -1
    * @param starts agents start locations - x, y coords
    * @param goals agents goal locations - x, y coords
    * @param scenstr name of the scenario, for logging purposes, default is "scen"
    * @param mapstr name of the map, for logging purposes, default is "map"
    */
    _MAPFSAT_Instance(std::vector<std::vector<int> >&, std::vector<std::pair<int,int> >&, std::vector<std::pair<int,int> >&, std::string = "scen", std::string = "map");

    /** Set the number of agents to be computed.
    *
    * Computes reachability and lower bounds for agents that were not used before. Does not check if ags is out of bounds.
    *
    * @param ags number of agents.
    */
    void SetAgents(int);

	/** Set locations to be avoided.
    *
    * Stores the provided argument.
    *
    * @param avoid x,y coords in time to be avoided by all agents.
    */
	void LoadAvoidData(std::vector<std::pair<std::pair<int,int>,int> >&);
	
	int GetMksLB(size_t);
	int GetSocLB(size_t);

	_MAPFSAT_Vertex IDtoCoords(int);
	bool HasNeighbor(_MAPFSAT_Vertex, int);
	bool HasNeighbor(int, int);
	int GetNeighbor(int, int);
	
	int FirstTimestep(int, int);
	int LastTimestep(int, int, int, int, int);
	int OppositeDir(int);

	void DebugPrint(std::vector<std::vector<int> >&);
	void DebugPrint(std::vector<int>&);
	void DebugPrint(std::vector<_MAPFSAT_Vertex>&);

	std::vector<_MAPFSAT_Agent> agents;
	std::vector<std::vector<int> > length_from_start;
	std::vector<std::vector<int> > length_from_goal;
	std::vector<int> SP_lengths;
	std::vector<_MAPFSAT_Avoid> avoid_locations;

	std::vector<std::vector<int> > map;
	size_t height;
	size_t width;
	size_t number_of_vertices;

	std::string scen_name;
	std::string map_name;

private:
	void LoadAgents(std::string, std::string);
	void LoadAgentsData(std::vector<std::pair<int,int> >&, std::vector<std::pair<int,int> >&);
	void LoadMap(std::string);
	void LoadMapData(std::vector<std::vector<int> >&);
	void BFS(std::vector<int>&, _MAPFSAT_Vertex);
	
	std::vector<int> mks_LBs;
	std::vector<int> soc_LBs;
	std::vector<_MAPFSAT_Vertex> coord_list;

	size_t last_number_of_agents;
};

#endif
#ifndef _logger_h_INCLUDED
#define _logger_h_INCLUDED


class _MAPFSAT_Logger
{

public:
	/** Constructor of _MAPFSAT_Logger.
    *
    * Stores all relevant data about the solution.
    *
    * @param instance pointer to a _MAPFSAT_Instance.
    * @param enc name of the encoding used by the soler.
    * @param log_level logging details - 0 = no log, 1 = inline log, 2 = human readable log. Default is 0.
    * @param log_file file output of the log. If "", print to stdout. Default is "".
    */
    _MAPFSAT_Logger(_MAPFSAT_Instance*, std::string, int = 0, std::string = "");

    /** Resets all currently stored data and prepares for the next solver call.
    *
    * @param ags number of agents to be solved in the next solver call.
    */
    void NewInstance(int);

    /** Prints stored statistics.
    *
    * Print is formatted as specified in the constructor. If there is no solution due to timeout, the log may not be correct.
    *
    */
    void PrintStatistics();

	std::string map_name;
	std::string scen_name;
	std::string encoding;
	int solution_mks;
	int solution_soc;
	int building_time;
	int solving_time;
	int nr_vars;
	int nr_clauses;
	int solver_calls;
	int agents;
	int mksLB;
	int socLB;	
	int res;
	
private:
	_MAPFSAT_Instance* inst;
	std::string log_file;

	int print_type;

};

#endif
#ifndef _solver_common_h_INCLUDED
#define _solver_common_h_INCLUDED

#include <unistd.h>
#include <chrono>
#include <thread>
#include <unordered_map>
#include <tuple>


struct _MAPFSAT_TEGAgent
{
	int first_variable;
	int first_timestep;
	int last_timestep;
};

struct _MAPFSAT_Shift
{
	int first_varaible;
	std::vector<int> timestep;
};

class _MAPFSAT_ISolver
{
public:
    virtual ~_MAPFSAT_ISolver() {};
    
    /** Perform a solve call.
    * 
    * Creates CNF formula based on the specified encoding and calls SAT solver to find a solution.
    *
    * @param ags number of agents in the solve call.
    * @param delta the initial delta. Default is 0.
    * @param oneshot option to perform just one solver call with the given delta without incrementing. Default is false.
	* @param keep_plan option save the found plan. The found plan can be retrieved by GetPlan function. Default is false.
    */
    int Solve(int, int = 0, bool = false, bool = false);

    /** Set data before solving.
    * 
    * Should be performed before the first solve. The stored data will be remembered for all of the solve calls.
    *
    * @param instance pointer to a _MAPFSAT_Instance.
    * @param logger pointer to a _MAPFSAT_Logger.
    * @param timeout timeout for each solve call in [s]. Gets reseted after each solve call.
	* @param CNF_file file name to print the created CNF formula. If no file is specified, the formula is not printed. Default is "".
    * @param quiet option to suppress any print to stdout. Default is false.
    * @param print_paths option to print found paths. Default is false.
	* @param use_avoids option to avoid certain postions in time. The avoid data is stored in the _MAPFSAT_Instance class. Default is false.
    */
    void SetData(_MAPFSAT_Instance*, _MAPFSAT_Logger*, int, std::string = "", bool = false, bool = false, bool = false);

	/** Returns the found plan.
    * 
    * Returns the found plan if the solve was success and keep_plan was set to true
    *
    */
	std::vector<std::vector<int> > GetPlan();

protected:
	_MAPFSAT_Instance* inst;
	_MAPFSAT_Logger* log;
	std::string solver_name;
	int cost_function; // 1 = mks, 2 = soc
	int movement; // 1 = parallel, 2 = pebble
	int lazy_const; // 1 = all at once, 2 = lazy
	int timeout;
	std::string cnf_file;
	bool quiet;
	bool print_plan;
	bool use_avoid;
	bool keep_plan;
	int solver_to_use = 1; // 1 = kissat, 2 = monosat

	int agents;
	int vertices;
	int delta;
	int max_timestep;
	std::stringstream cnf_printable;

	_MAPFSAT_TEGAgent** at;
	_MAPFSAT_TEGAgent*** pass;
	_MAPFSAT_Shift** shift;
	int* shift_times_start;
	int* shift_times_end;
	int at_vars;

	int nr_vars;
	int nr_clauses;
	int solver_calls;

	std::vector<std::vector<int> > plan;

	bool conflicts_present;

	std::vector<std::tuple<int,int,int,int> > vertex_conflicts;
	std::vector<std::tuple<int,int,int,int,int> > swap_conflicts;
	std::vector<std::tuple<int,int,int,int,int> > pebble_conflicts;

	// solver
	void* SAT_solver;

	// before solving
	void PrintSolveDetails(int);

	// virtual encoding to be used
	virtual int CreateFormula(int) = 0;

	// creating formula
	int CreateAt(int, int);
	int CreatePass(int, int);
	int CreateShift(int, int);

	void CreatePossition_Start();
	void CreatePossition_Goal();
	void CreatePossition_NoneAtGoal();

	void CreateConf_Vertex();
	void CreateConf_Swapping_At();
	void CreateConf_Swapping_Pass();
	void CreateConf_Swapping_Shift();
	void CreateConf_Pebble_At();
	void CreateConf_Pebble_Pass();
	void CreateConf_Pebble_Shift();

	void CreateConf_Vertex_OnDemand();
	void CreateConf_Swapping_At_OnDemand();
	void CreateConf_Swapping_Pass_OnDemand();
	void CreateConf_Swapping_Shift_OnDemand();
	void CreateConf_Pebble_At_OnDemand();
	void CreateConf_Pebble_Pass_OnDemand();
	void CreateConf_Pebble_Shift_OnDemand();

	void CreateMove_NoDuplicates();
	void CreateMove_NextVertex_At();
	void CreateMove_EnterVertex_Pass();
	void CreateMove_LeaveVertex_Pass();
	void CreateMove_ExactlyOne_Shift();
	void CreateMove_NextVertex_Shift();

	void CreateMove_Graph_MonosatPass();

	int CreateConst_LimitSoc(int);
	void CreateConst_Avoid();

	// solving
	int InvokeSolver_Kissat(int);
	int InvokeSolver_Monosat(int);
	
	// auxiliary SAT solver functions
	void CreateSolver();
	int InvokeSolver(int);
	static void WaitForTerminate(int, void*, bool&);
	int VarToID(int, bool, int&, std::unordered_map<int, int>&);
	void AddClause(std::vector<int>);

	// plan outputting functions
	int NormalizePlan();
	void VerifyPlan();
	void GenerateConflicts();

	// cleanup functions
	bool TimesUp(std::chrono::time_point<std::chrono::high_resolution_clock>, std::chrono::time_point<std::chrono::high_resolution_clock>, int);
	void CleanUp(bool);
};

/******************************************************************/
/*********************** Specific Encodings ***********************/
/******************************************************************/

class _MAPFSAT_AtParallelMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtParallelMksAll(std::string name = "at_parallel_mks_all");
	~_MAPFSAT_AtParallelMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtParallelSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtParallelSocAll(std::string name = "at_parallel_soc_all");
	~_MAPFSAT_AtParallelSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtPebbleMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtPebbleMksAll(std::string name = "at_pebble_mks_all");
	~_MAPFSAT_AtPebbleMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtPebbleSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtPebbleSocAll(std::string name = "at_pebble_soc_all");
	~_MAPFSAT_AtPebbleSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassParallelMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassParallelMksAll(std::string name = "pass_parallel_mks_all");
	~_MAPFSAT_PassParallelMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassParallelSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassParallelSocAll(std::string name = "pass_parallel_soc_all");
	~_MAPFSAT_PassParallelSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassPebbleMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassPebbleMksAll(std::string name = "pass_pebble_mks_all");
	~_MAPFSAT_PassPebbleMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassPebbleSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassPebbleSocAll(std::string name = "pass_pebble_soc_all");
	~_MAPFSAT_PassPebbleSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftParallelMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftParallelMksAll(std::string name = "shift_parallel_mks_all");
	~_MAPFSAT_ShiftParallelMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftParallelSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftParallelSocAll(std::string name = "shift_parallel_soc_all");
	~_MAPFSAT_ShiftParallelSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftPebbleMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftPebbleMksAll(std::string name = "shift_pebble_mks_all");
	~_MAPFSAT_ShiftPebbleMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftPebbleSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftPebbleSocAll(std::string name = "shift_pebble_soc_all");
	~_MAPFSAT_ShiftPebbleSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtParallelMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtParallelMksLazy(std::string name = "at_parallel_mks_lazy");
	~_MAPFSAT_AtParallelMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtParallelSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtParallelSocLazy(std::string name = "at_parallel_soc_lazy");
	~_MAPFSAT_AtParallelSocLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtPebbleMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtPebbleMksLazy(std::string name = "at_pebble_mks_lazy");
	~_MAPFSAT_AtPebbleMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_AtPebbleSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_AtPebbleSocLazy(std::string name = "at_pebble_soc_lazy");
	~_MAPFSAT_AtPebbleSocLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassParallelMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassParallelMksLazy(std::string name = "pass_parallel_mks_lazy");
	~_MAPFSAT_PassParallelMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassParallelSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassParallelSocLazy(std::string name = "pass_parallel_soc_lazy");
	~_MAPFSAT_PassParallelSocLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassPebbleMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassPebbleMksLazy(std::string name = "pass_pebble_mks_lazy");
	~_MAPFSAT_PassPebbleMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_PassPebbleSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_PassPebbleSocLazy(std::string name = "pass_pebble_soc_lazy");
	~_MAPFSAT_PassPebbleSocLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftParallelMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftParallelMksLazy(std::string name = "shift_parallel_mks_lazy");
	~_MAPFSAT_ShiftParallelMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftParallelSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftParallelSocLazy(std::string name = "shift_parallel_soc_lazy");
	~_MAPFSAT_ShiftParallelSocLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftPebbleMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftPebbleMksLazy(std::string name = "shift_pebble_mks_lazy");
	~_MAPFSAT_ShiftPebbleMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_ShiftPebbleSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_ShiftPebbleSocLazy(std::string name = "shift_pebble_soc_lazy");
	~_MAPFSAT_ShiftPebbleSocLazy() {};
private:
	int CreateFormula(int);
};

/*****************************************************************/
/*********************** Monosat Encodings ***********************/
/*****************************************************************/

class _MAPFSAT_MonosatParallelMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatParallelMksAll(std::string name = "monosat_parallel_mks_all");
	~_MAPFSAT_MonosatParallelMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatParallelSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatParallelSocAll(std::string name = "monosat_parallel_soc_all");
	~_MAPFSAT_MonosatParallelSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatPebbleMksAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatPebbleMksAll(std::string name = "monosat_pebble_mks_all");
	~_MAPFSAT_MonosatPebbleMksAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatPebbleSocAll : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatPebbleSocAll(std::string name = "monosat_pebble_soc_all");
	~_MAPFSAT_MonosatPebbleSocAll() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatParallelMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatParallelMksLazy(std::string name = "monosat_parallel_mks_lazy");
	~_MAPFSAT_MonosatParallelMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatParallelSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatParallelSocLazy(std::string name = "monosat_parallel_soc_lazy");
	~_MAPFSAT_MonosatParallelSocLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatPebbleMksLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatPebbleMksLazy(std::string name = "monosat_pebble_mks_lazy");
	~_MAPFSAT_MonosatPebbleMksLazy() {};
private:
	int CreateFormula(int);
};

class _MAPFSAT_MonosatPebbleSocLazy : public _MAPFSAT_ISolver
{
public:
	_MAPFSAT_MonosatPebbleSocLazy(std::string name = "monosat_pebble_soc_lazy");
	~_MAPFSAT_MonosatPebbleSocLazy() {};
private:
	int CreateFormula(int);
};

#endif
