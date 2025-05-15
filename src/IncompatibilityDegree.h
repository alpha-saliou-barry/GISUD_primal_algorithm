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
	std::map<int, std::string> nodes_columns;
	// Vector of tasks corresponding to task order
	std::vector<std::string> tasks_order_;

public:
	// Get the node number identified by "task_" and column "column"
	int getNode(std::string task, std::string column);
	

	// Constructor of incompatibility degree
	// "problem" is a pointer on the problem
	// "positive_columns" is the vector of columns in P
	// "tasks_order" is the tasks order
	IncompatibilityDegree(ISUD_Base* problem, std::vector <IB_Column*> positive_columns, std::vector < std::string > tasks);
	
	// Empty the incompatibility degree graph
	void emptyGraph();
	

	// Construct the graph to compute inc degree approximation
	// "column" is the column for which we want to construct the graph
	// "contribsOut" is a map of additional contributions, can be NULL
	bool constructGraph(IB_Column* column, std::unordered_map < std::string, int >* contribsOut = NULL);
	
	// Return the incompatibility degree of column "column"
	// "constribsOut" is additional contributions
	// "involvedColumns" is the set of columns indices making the column "column" binary compatible
	int getIncompatibilityDegree(IB_Column*, std::unordered_map < std::string, int >* contribsOut = NULL, std::set<std::string>* involvedColumns = NULL);

	// Compute HASH for a path in the graph
	std::string computePathString(std::pair<std::vector<int>, int> path);
};
