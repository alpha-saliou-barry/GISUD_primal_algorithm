#include "IncompatibilityDegree.h"
#include <string>   

// Constructor of incompatibility degree
IncompatibilityDegree::IncompatibilityDegree(ISUD_Base* problem, std::vector <IB_Column*> positive_columns,
	std::vector <std::string> tasks_order, std::map<std::pair<std::string, std::string>, std::pair<int, int>>* tasksCorrespondance_): tasks_order_(tasks_order) {
	
	if(tasksCorrespondance_ != NULL) {
		tasksCorrespondance = *tasksCorrespondance_;
	}
	// Initialise le calculateur de degré d'incompatibilité
	// Prend en paramètre les colonnes positives ainsi que l'ordre sur les tâches
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

// Get sets of nb_columns in positive_columns
std::vector < std::set < std::string > > IncompatibilityDegree::getSets(std::vector < std::string > positive_columns, int nb_columns) {
	// Retourne tous les ensembles de nb_columns dans l'ensemble des positive_columns
	if (nb_columns > positive_columns.size()) {
		return std::vector < std::set < std::string > >();
	}
	std::vector < std::set < std::string > > parts;
	if (nb_columns == 1) {
		for (auto element : positive_columns) {
			std::set < std::string > part = std::set < std::string >();
			part.insert(element);

			parts.push_back(part);
		}
		
		return parts;
	}



	for (auto value_ : positive_columns) {
		std::vector < std::string > new_set;
		for (auto val_ : positive_columns) {
			if (val_ > value_) {
				new_set.push_back(val_);
			}
		}

		for (auto part : getSets(new_set, nb_columns - 1)) {
			part.insert(value_);
			parts.push_back(part);
		}

	}

	return parts;
}

// Get the node number identified by task_ and columns as columns set
int IncompatibilityDegree::getNode(std::string task_, std::set < std::string > columns) {
	// Retourne le numéro du noeud identifié par task_ et columns comme ensemble de colonnes
	std::string second_string;
	std::string first_col;
	if (columns.empty()) {
		second_string = "NIL";
	}
	else {
		std::vector < std::string > other_parts;
		for (auto column : columns) {
			other_parts.push_back(column);
			first_col = column;
		}

		std::sort(other_parts.begin(), other_parts.end());
		for (auto part : other_parts) {
			second_string += part + "-";
		}
	}

	std::string node_str = task_ + "|" + second_string;
	if (nodes.find(node_str) == nodes.end()) {
		if(tasksCorrespondance.size()) {
			nodes_tasks[nb_nodes] = tasksCorrespondance[std::pair<std::string, std::string>(task_, first_col)];
		}
		nodes[node_str] = nb_nodes;
		nodes_columns[nb_nodes] = columns;
		
		nb_nodes += 1;
	}


	return nodes[node_str];
}

// Construct the graph to compute inc degree approximation
bool IncompatibilityDegree::constructGraph(IB_Column* column, std::unordered_map < std::string, int >* contribsOut) {
	// Construit le graphe permettant de calculer le degré d'incompatibilité
	
	// On construit les noeuds du graphes
	getNode("S", std::set<std::string>());
	getNode("T", std::set<std::string>());

	// On construit les arêtes du graphe
	std::vector < std::pair < std::string, std::set <std::string> > > previous_nodes;
	previous_nodes.push_back(std::pair < std::string, std::set < std::string > >("S", std::set < std::string >()));
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
		std::cout << "pb2 " << std::endl;
		return false;
	}

	tasks_pairs.push_back("T");
	for (auto task_pair : tasks_pairs) {
		std::vector < std::pair < std::string, std::set < std::string> > > new_nodes;
		if (task_pair == "T") {
			new_nodes.push_back(std::pair < std::string, std::set < std::string > >("T", std::set < std::string >()));
		}
		else {
			for (auto set_ : getSets(positives_columns_per_task[task_pair], 1)) {
				new_nodes.push_back(std::pair < std::string, std::set < std::string > >(task_pair, set_));
			}
		}


		for (auto previous_node : previous_nodes) {
			for (auto new_node : new_nodes) {
				int d = 0;
				if (new_node.first != "T") {
					for (auto new_column : new_node.second) {
						if ((predecessors[new_column][new_node.first] != "NIL") && (previous_node.second.find(new_column) == previous_node.second.end()
							|| predecessors[new_column][new_node.first] != previous_node.first)) {
							d += 1;
						}
					}
				}

				if (previous_node.first != "S") {
					for (auto previous_column : previous_node.second) {
						if ((successors[previous_column][previous_node.first] != "NIL") && (new_node.second.find(previous_column) == new_node.second.end()
							|| successors[previous_column][previous_node.first] != new_node.first)) {
							d += 1;
						}
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


// Return the incompatibility degree of columns column
int IncompatibilityDegree::getIncompatibilityDegree(IB_Column* column, std::unordered_map < std::string, int >* contribsOut, std::set<std::string>* involvedColumns) {
	// Calcule le degré d'incompatibilité de la colonne column dans le graphe. ContribsOut est vide, involvedColumns est l'adresse à laquelle vont être stockée les colonnes 
	// rendant compatible binaire une colonne de degré d'incompatibilité 0
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
				for (auto column : nodes_columns[predecessor]) {
					involvedColumns->insert(column);
				}
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

// Return all paths of the column lower than phase phase
std::vector<std::pair<std::map<int, int>, int>> IncompatibilityDegree::getAllPaths(IB_Column* column, int phase) {
	emptyGraph();
	constructGraph(column, NULL);
	
	/*
	for (auto node_pair : nodes) {
		int node = node_pair.second;
		for (auto edge : edges[node]) {
			std::cout << "Edge : " << node << " > " << edge.first  << " : " << edge.second << std::endl;
		}
	}*/

	// Enumeration de tous les chemins
	std::map<int, std::vector<std::pair<std::vector<int>, int>>> nodesAllPaths;
	std::map<int, std::map<std::string, int>> nodesPathsStrings;

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
	int puit = nodes["T|NIL"];
	nodesAllPaths[nodes["S|NIL"]] =  std::vector<std::pair<std::vector<int>, int>>({ std::pair<std::vector<int>, int>({nodes["S|NIL"]}, 0) });
	bool finished = false;
	std::map<int, int> predecessors_djikstra;

	//std::cout << 0 << " > " << puit << std::endl;

	while (!finished) {
		std::vector < int > c_layer;
		for (auto node : current_layer) {
			std::vector<std::pair<std::vector<int>, int>> nodePaths = nodesAllPaths[node];
			std::vector < std::pair < int, int > > next_nodes = edges[node];

			for (auto next_node : next_nodes) {
				if (nodesAllPaths.find(next_node.first) == nodesAllPaths.end()) {
					nodesAllPaths[next_node.first] = std::vector<std::pair<std::vector<int>, int>>();
					nodesPathsStrings[next_node.first] = std::map<std::string, int>();
				}

				for (auto path : nodePaths) {
					std::vector<int> pathContent = path.first;
					int pathDistance = path.second;
					pathDistance += next_node.second;
					if (pathDistance <= phase) {
						std::string pathString = computePathString(path);
						pathContent.push_back(next_node.first);
						if (nodesPathsStrings[next_node.first].find(pathString) != nodesPathsStrings[next_node.first].end()) {
							int path_id = nodesPathsStrings[next_node.first][computePathString(path)];

							if (nodesAllPaths[next_node.first][path_id].second > pathDistance) {
								nodesAllPaths[next_node.first][path_id] = std::pair<std::vector<int>, int>(pathContent, pathDistance);
							}
						}
						else {
							nodesAllPaths[next_node.first].push_back(std::pair<std::vector<int>, int>(pathContent, pathDistance));
							nodesPathsStrings[next_node.first][pathString] = nodesAllPaths[next_node.first].size()  - 1;
						}
					}
					


				}

				numberOfEdges[next_node.first] -= 1;
				if (numberOfEdges[next_node.first] == 0) {
					//std::cout << "Fini : " << next_node.first << std::endl;
					c_layer.push_back(next_node.first);
					if (next_node.first == puit) {
						finished = true;
					}
				}
			}
		}

		current_layer = c_layer;
	}

	/*
	for (auto pair : nodesAllPaths) {
		std::cout << pair.first << " | " << pair.second.size() << std::endl;
	}

		for (auto path : nodesAllPaths[puit]) {
			std::cout << "Nouvelle colonne : " << std::endl;
			for (auto node : path.first) {
				if (node != 0 && node != puit) {
					std::cout << nodes_infos[node].first << std::endl;
				}
			}
			std::cout << "Degre d incompatibilite : " << path.second << std::endl;
		}
		*/
	std::vector<std::pair<std::map<int, int>, int>> returnValue;
	for (auto path : nodesAllPaths[puit]) {
		std::map<int, int> columnContribs;
		for (auto node : path.first) {
			if (node != nodes["S|NIL"] && node != puit) {
				columnContribs[nodes_tasks[node].first] = nodes_tasks[node].second;
				//std::cout << nodes_tasks[node].first << " " << nodes_tasks[node].second << " | ";
			}
		}

		//std::cout << std::endl;



		returnValue.push_back(std::pair<std::map<int, int>, int>(columnContribs, path.second));
	}

	//std::cout << "t((((((((((((((((((((((((((()))))))))))))))))))))))))))" << std::endl;

	if(returnValue.size() == 0) {
		std::cout << "bizarre" << std::endl;
	}

	return returnValue;
}