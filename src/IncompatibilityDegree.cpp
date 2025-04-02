#include "IncompatibilityDegree.h"
#include <string>   

// Constructor, "problem" is the problem, "positive_columns" is the columns that are in the current_solution
// "tasks" is the ordered tasks

IncompatibilityDegree::IncompatibilityDegree(ISUD_Base* problem, std::vector <IB_Column*> positive_columns,
	std::vector <std::string> tasks_order): tasks_order_(tasks_order) {
	psolutionMethod_ = problem;
	int i = 0;
	std::map<std::string, int> num;
	for (auto  task : tasks_order) {
		positives_columns_per_task[task] = std::vector < std::string >();
		num[task] = i;
		i += 1;
	}


	for (auto column : positive_columns) {
		std::string current_predecessor = "NIL";
		std::map <std::string, std::string> column_predecessors;
		std::map <std::string, std::string> column_sucessors;

		for (auto task : column->getContribsOrder()) {
			if (column->findContribution(task) != 0) {
				
				if (rhs.find(task) == rhs.end()) {
					rhs[task] = 0;
				}

				rhs[task] += 1;
				positives_columns_per_task[task].push_back(column->getName());

				column_predecessors[task] = current_predecessor;
				if (current_predecessor != "NIL") {
					column_sucessors[current_predecessor] = task;
				}

				current_predecessor = task;
			}
		}

		


		if (current_predecessor != "NIL") {
			column_sucessors[current_predecessor] = "NIL";
		}

		predecessors[column->getName()] = column_predecessors;
		successors[column->getName()] = column_sucessors;
	}
}

// Clear the incompatibility degree graph
void IncompatibilityDegree::emptyGraph() {
	nb_nodes = 0;
	nodes = std::map < std::string, int >();
	edges = std::map < int, std::vector < std::pair < int, int> > >();
}

// Return the node identifier associated to task "task" and column "column"
int IncompatibilityDegree::getNode(std::string task_, std::string column) {
	if (column == "") {
		column = "NIL";
	}

	std::string node_str = task_ + "|" + column;
	if (nodes.find(node_str) == nodes.end()) {
		nodes[node_str] = nb_nodes;
		if (column != "NIL") {
			nodes_column[nb_nodes] = column;
		}
		
		nb_nodes += 1;
	}


	return nodes[node_str];
}

// Construct the incompatibility degree graph for column "column"
bool IncompatibilityDegree::constructGraph(IB_Column* column, std::unordered_map < std::string, int >* contribsOut) {
	// Nodes construction
	getNode("S", "NIL");
	getNode("T", "NIL");

	// Edges construction
	std::vector < std::pair < std::string, std::string> > previous_nodes;
	previous_nodes.push_back(std::pair < std::string, std::string>("S", "NIL"));
	std::vector < std::string> ordered_tasks;
	std::map <std::string, int > contributions;
	for (auto task : column->getContribsOrder()) {
		if (column->findContribution(task) != 0 && rhs[task] > 0) {
			ordered_tasks.push_back(task);
			contributions[task] = 1;
		}
		else if (column->findContribution(task) != 0 && rhs[task] == 0) {
			return false;
		}
	}

	if (contribsOut != NULL) {
		for (auto bnz : *contribsOut) {
			if (std::find(ordered_tasks.begin(), ordered_tasks.end(), bnz.first) == ordered_tasks.end()) {
				ordered_tasks.push_back(bnz.first);
				contributions[bnz.first] = bnz.second;
			}
			else {
				contributions[bnz.first] += bnz.second;
			}
		}
	}

	for (auto contribs_pair : contributions) {
		if (contribs_pair.second > rhs[contribs_pair.first]) {
			return false;
		}

		if (contribs_pair.second <= 0) {
			std::vector<std::string>::iterator position = std::find(ordered_tasks.begin(), ordered_tasks.end(), contribs_pair.first);
			if (position != ordered_tasks.end()) {
				ordered_tasks.erase(position);
			}
		}
	}

	if (ordered_tasks.size() == 0) {
		return false;
	}

	ordered_tasks.push_back("T");
	for (auto task_pair : ordered_tasks) {
		std::vector < std::pair < std::string, std::string > > new_nodes;
		if (task_pair == "T") {
			new_nodes.push_back(std::pair < std::string, std::string >("T", "NIL"));
		}
		else {
			for (auto column : positives_columns_per_task[task_pair]) {
				new_nodes.push_back(std::pair < std::string, std::string >(task_pair, column));
			}
		}


		for (auto previous_node : previous_nodes) {
			for (auto new_node : new_nodes) {
				int d = 0;
				if (new_node.first != "T") {
					std::string new_column = new_node.second;
					
					if ((predecessors[new_column][new_node.first] != "NIL") && (previous_node.second != new_column
							|| predecessors[new_column][new_node.first] != previous_node.first)) {
							d += 1;
						}
					}

				if (previous_node.first != "S") {

						std::string previous_column = previous_node.second;
						if ((successors[previous_column][previous_node.first] != "NIL") && (new_node.second != previous_column
							|| successors[previous_column][previous_node.first] != new_node.first)) {
							d += 1;
						}
					}

				int first_node = getNode(previous_node.first, previous_node.second);
				int second_node = getNode(new_node.first, new_node.second);
				if (edges.find(first_node) == edges.end()) {
					edges[first_node] = std::vector < std::pair< int, int > >();
				}

				edges[first_node].push_back(std::pair < int, int >(second_node, d));
			}
		}

		previous_nodes = new_nodes;
	}

	return true;

}

// Compute the incompatibility degree of column "column"
// If column is binary compatible, return the positive columns that make "column" compatible in involvedColumns
int IncompatibilityDegree::getIncompatibilityDegree(IB_Column* column, std::unordered_map < std::string, int >* contribsOut, std::set<std::string>* involvedColumns) {
	emptyGraph();
	bool graphConstructed = constructGraph(column, contribsOut);
	if (!graphConstructed) {
		return 6;
	}

	std::map < int, int > numberOfEdges;
	std::map < int, int > distances;

	for (auto edge_pair : edges) {
		for (auto edge_pair_2 : edge_pair.second) {
			numberOfEdges[edge_pair_2.first] += 1;
		}
	}

	std::vector < int > current_layer;
	current_layer.push_back(nodes["S|NIL"]);
	distances[nodes["S|NIL"]] = 0;
	bool finished = false;
	std::map<int, int> predecessors_djikstra;
	while (!finished) {
		std::vector < int > c_layer;
		for (auto node : current_layer) {
			std::vector < std::pair < int, int > > next_nodes = edges[node];

			for (auto next_node : next_nodes) {
				if (distances.find(next_node.first) == distances.end() || distances[next_node.first] > distances[node] + next_node.second) {
					distances[next_node.first] = distances[node] + next_node.second;
					predecessors_djikstra[next_node.first] = node;
				}

				numberOfEdges[next_node.first] -= 1;
				if (numberOfEdges[next_node.first] == 0) {
					c_layer.push_back(next_node.first);
					if (next_node.first == nodes["T|NIL"]) {
						finished = true;
					}
				}
			}
		}

		current_layer = c_layer;
	}

	if (involvedColumns != NULL) {
		int current_node = nodes["T|NIL"];
		int final_node = nodes["S|NIL"];
		while (current_node != final_node) {
			int predecessor = predecessors_djikstra[current_node];

			if (predecessor != nodes["S|NIL"] && predecessor != nodes["T|NIL"]) {
				involvedColumns->insert(nodes_column[predecessor]);
			}
			

			current_node = predecessor;
		}
	}
	

	return distances[nodes["T|NIL"]];
}
