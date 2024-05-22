#include <vector>
#include <unordered_map>
#include <algorithm>
#include "ISUD_Base.h"
#include "IB_Column.h"
#include <map> 
#include <set>

class IncompatibilityDegree {
protected:
	// Task name to task order
	std::map<std::string, int> tasksOrderInv;
	// Tasks correspondance
	std::map<std::pair<std::string, std::string>, std::pair<int, int>> tasksCorrespondance;
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;
	// RHS of the tasks
	std::unordered_map < std::string, int > rhs;
	// Nodes hash table 
	std::map < std::string, int > nodes;
	// Edges map
	std::map < int, std::vector < std::pair < int, int > > > edges;
	// Total number of node
	int nb_nodes;
	// Predecessors of task 
	std::map < std::string, std::map <std::string, std::string> > predecessors;
	// Successors of task
	std::map < std::string, std::map <std::string, std::string> > successors;
	// Positive columns per task
	std::map < std::string, std::vector < std::string > > positives_columns_per_task;
	// Map node number to task
	std::map<int, std::pair<int, int>> nodes_tasks;
	// Map nodes number to columns covered by node
	std::map<int, std::set<std::string>> nodes_columns;
	// Disaggregated tasks no
	std::set<int> disaggregatedTasks;
	// Correspondance between initial task and disaggregated task
	std::map<int, int> correspondance;
	// Vector of tasks corresponding to task order
	std::vector<std::string> tasks_order_;

public:
	// Get sets of nb_columns in positive_columns
	std::vector < std::set < std::string > > getSets(std::vector < std::string > parts, int k);
	
	// Get the node number identified by task_ and columns as columns set
	int getNode(std::string task, std::set < std::string > pc);
	
	// Constructor of incompatibility degree
	IncompatibilityDegree(ISUD_Base* problem, std::vector <IB_Column*> positive_columns, std::vector < std::string > tasks,
	std::map<std::pair<std::string, std::string>, std::pair<int, int>>* tasksCorrespondance_ = NULL);
	
	// Empty the incompatibility degree graph
	void emptyGraph();
	
	// Construct the graph to compute inc degree approximation
	bool constructGraph(IB_Column* column, std::unordered_map < std::string, int >* contribsOut = NULL);
	
	// Return the incompatibility degree of columns column
	int getIncompatibilityDegree(IB_Column*, std::unordered_map < std::string, int >* contribsOut = NULL, std::set<std::string>* involvedColumns = NULL);
	
	// Return all paths of the column lower than phase phase
	std::vector<std::pair<std::map<int, int>, int>> getAllPaths(IB_Column* column, int phase);

	// Compute HASH for a path in the graph
	std::string computePathString(std::pair<std::vector<int>, int> path);
};
