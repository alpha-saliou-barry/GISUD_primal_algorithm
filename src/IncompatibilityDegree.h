#include <vector>
#include <unordered_map>
#include <algorithm>
#include "ISUD_Base.h"
#include "IB_Column.h"
#include <map> 
#include <set>

class IncompatibilityDegree {
protected:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;

	// For each tasks, gives the rhs (b_i)
	std::unordered_map < std::string, int > rhs;

	// String node to int
	std::map < std::string, int > nodes;

	// Edges of graph
	std::map < int, std::vector < std::pair < int, int > > > edges;

	// Number of nodes
	int nb_nodes;

	// Predecessors in column of covered tasks
	std::map < std::string, std::map <std::string, std::string> > predecessors;

	// Successors in column of covered tasks
	std::map < std::string, std::map <std::string, std::string> > successors;

	// For each task, it gives the positive columns that cover the task
	std::map < std::string, std::vector < std::string > > positives_columns_per_task;

	// For each node, gives the column associated to node
	std::map<int, std::string> nodes_column;

	// Tasks order
	std::vector<std::string> tasks_order_;

public:
	// Return the node identifier associated to task "task" and column "column"
	int getNode(std::string task, std::string column);

	// Constructor, "problem" is the problem, "positive_columns" is the columns that are in the current_solution
	// "tasks" is the ordered tasks

	IncompatibilityDegree(ISUD_Base* problem, std::vector <IB_Column*> positive_columns, std::vector < std::string > tasks);
	// Clear the incompatibility degree graph
	void emptyGraph();
	// Construct the incompatibility degree graph for column "column"
	bool constructGraph(IB_Column* column, std::unordered_map < std::string, int >* contribsOut = NULL);
	// Compute the incompatibility degree of column "column"
	// If column is binary compatible, return the positive columns that make "column" compatible in involvedColumns
	int getIncompatibilityDegree(IB_Column* column, std::unordered_map < std::string, int >* contribsOut = NULL, std::set<std::string>* involvedColumns = NULL);
	
};
