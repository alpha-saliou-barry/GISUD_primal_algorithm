#include "IncompatibilityDegree.h"
#include <string>   

// Constructor of incompatibility degree
// "problem" is a pointer on the problem
// "positive_columns" is the vector of columns in P
// "tasks_order" is the tasks order
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

// Empty the incompatibility degree graph
void IncompatibilityDegree::emptyGraph() {
	// Vide la structure de graphe permettant de calculer le degré d'incompatibilité
	nb_nodes = 0;
	nodes = std::map < std::string, int >();
	edges = std::map < int, std::vector < std::pair < int, int> > >();
}

// Get the node number identified by "task_" and column "column"
int IncompatibilityDegree::getNode(std::string task_, std::string column) {
	std::string second_string;
	std::string first_col;
	if (column == "") {
		second_string = "NIL";
	}

	else {
		second_string += column;
	}

	std::string node_str = task_ + "|" + second_string;
	if (nodes.find(node_str) == nodes.end()) {
		nodes[node_str] = nb_nodes;
		nodes_columns[nb_nodes] = column;
		
		nb_nodes += 1;
	}


	return nodes[node_str];
}

// Construct the graph to compute inc degree approximation
// "column" is the column for which we want to construct the graph
// "contribsOut" is a map of additional contributions, can be NULL
bool IncompatibilityDegree::constructGraph(IB_Column* column, std::unordered_map < std::string, int >* contribsOut) {
	// Construit le graphe permettant de calculer le degré d'incompatibilité
	
	// On construit les noeuds du graphes
	getNode("S", "");
	getNode("T", "");

	// On construit les arêtes du graphe
	std::vector < std::pair < std::string, std::string > > previous_nodes;
	previous_nodes.push_back(std::pair < std::string, std::string>("S", ""));
	std::vector < std::string> tasks_pairs;
	std::map <std::string, int > contributions;
	for (auto task : column->getContribsOrder()) {
		if (column->findContribution(task) != 0 && rhs[task] > 0) {
			tasks_pairs.push_back(task);
			contributions[task] = 1;
		}
		else if (column->findContribution(task) != 0 && rhs[task] == 0) {
			return false;
		}
	}

	if (contribsOut != NULL) {
		for (auto bnz : *contribsOut) {
			if (std::find(tasks_pairs.begin(), tasks_pairs.end(), bnz.first) == tasks_pairs.end()) {
				tasks_pairs.push_back(bnz.first);
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
			std::vector<std::string>::iterator position = std::find(tasks_pairs.begin(), tasks_pairs.end(), contribs_pair.first);
			if (position != tasks_pairs.end()) {
				tasks_pairs.erase(position);
			}
		}
	}

	if (tasks_pairs.size() == 0) {
		return false;
	}

	tasks_pairs.push_back("T");

	//std::cout << "Before tasks pairs" << std::endl;
	for (auto task_pair : tasks_pairs) {
		std::vector < std::pair < std::string, std::string > > new_nodes;
		if (task_pair == "T") {
			new_nodes.push_back(std::pair < std::string, std::string >("T", ""));
		}
		else {
			for (auto column: positives_columns_per_task[task_pair]) {
				new_nodes.push_back(std::pair < std::string, std::string >(task_pair, column));
			}
		}


		for (auto previous_node : previous_nodes) {
			for (auto new_node : new_nodes) {
				int d = 0;
				//std::cout << "New node : " << new_node.first << ", " << new_node.second;
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

				//std::cout << "Length put on arcs" << std::endl;

				int first_node = getNode(previous_node.first, previous_node.second);
				int second_node = getNode(new_node.first, new_node.second);
				if (edges.find(first_node) == edges.end()) {
					edges[first_node] = std::vector < std::pair< int, int > >();
				}

				edges[first_node].push_back(std::pair < int, int >(second_node, d));
				//std::cout << "Edges added" << std::endl;
			}
		}

		previous_nodes = new_nodes;
	}
	//std::cout << "After tasks pairs." << std::endl;

	return true;

}


// Return the incompatibility degree of column "column"
// "constribsOut" is additional contributions
// "involvedColumns" is the set of columns indices making the column "column" binary compatible
int IncompatibilityDegree::getIncompatibilityDegree(IB_Column* column, std::unordered_map < std::string, int >* contribsOut, std::set<std::string>* involvedColumns) {
	emptyGraph();
	//std::cout << "Graph construction" << std::endl;
	bool graphConstructed = constructGraph(column, contribsOut);
	//std::cout << "Graph constructed" << std::endl;
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
				involvedColumns->insert(nodes_columns[predecessor]);
			}
			

			current_node = predecessor;
		}
	}
	

	return distances[nodes["T|NIL"]];
}

// Compute HASH for a path in the graph
std::string IncompatibilityDegree::computePathString(std::pair<std::vector<int>, int> path) {
	std::string pathString;
	for (auto node : path.first) {
		if (node != nodes["S|NIL"] && node != nodes["T|NIL"]) {
			pathString += std::to_string(nodes_tasks[node].first) + ", " + std::to_string(nodes_tasks[node].second) + "|";
		}
	}
	return pathString;
}
