#include "ProblemReader.h"
#include <iostream>
#include <string>
#include <fstream>

// Read tasks from file
std::vector<std::pair<std::string, int>> read_tasks(std::string path) {
	std::ifstream file;
	file.open(path);

	bool start = false;
	bool in_task = false;
	std::string task_name = "";
	std::vector<std::pair<std::string, int>> return_result;
	while (file.good()) {
		std::string token = "";
		file >> token;
		if (token == "{") {
			start = true;
			in_task = false;
		}
		else if (token == "}") {
			start = false;
		}
		else if (start) {
			if (in_task) {
				return_result.push_back(std::pair<std::string, int>(task_name, std::stoi(token.substr(0, token.size() - 1))));
				in_task = false;
			}
			else {
				in_task = true;
				task_name = token;
			}
		}
	}

	return return_result;
}

// Read columns from path
std::vector<IB_Column*> read_columns(std::string path, double proba) {
	std::ifstream file;
	std::cout << path << std::endl;
	file.open(path);
	std::vector<IB_Column*> columns;
	while (file.good()) {
		std::string colName, intermediateToken, cost;
		file >> colName;
		if (!colName.size()) {
			continue;
		}
		file >> intermediateToken;
		file >> cost;
		file >> intermediateToken;
		std::string end_character = "";
		std::vector<std::pair<std::string, int>> contribs;
		while (end_character != ";") {
			std::string taskName, contrib;
			file >> taskName;
			file >> contrib;
			if (contrib.size() > 1) {
				end_character = contrib.substr(contrib.size() - 1, 1);
				contrib = contrib[0];
			}

			contribs.push_back(std::pair<std::string, int>(taskName, stoi(contrib)));
		}

		if (proba < 0 || true ) {
			columns.push_back(new IB_Column(colName, contribs));
			columns[columns.size() - 1]->setCost(std::stod(cost));
		}

	}

	return columns;
}

// Read initial solution from path
std::set<std::string> read_initial_solution(std::string path) {
	std::ifstream file;
	file.open(path);

	bool start = false;
	bool in_task = false;
	std::string task_name = "";
	std::set <std::string> return_result;
	while (file.good()) {
		std::string token = "";
		file >> token;
		if (token == "{") {
			start = true;
			in_task = false;
		}
		else if (token == "}") {
			start = false;
		}
		else if (start) {
			if (in_task) {
				in_task = false;
			}
			else {
				in_task = true;
				task_name = token;
				return_result.insert(task_name);
			}
		}
	}

	return return_result;
}


// Get random problem
ISUD_Base get_problem_random(std::string folder) {
	std::vector<std::pair<std::string, int>> tasks = read_tasks(folder + "/rhs.in");
	std::cout << "tasks read" << std::endl;
	std::vector<std::string> tasksNames;
	std::vector<int> rhs;
	for (auto pair : tasks) {
		tasksNames.push_back(pair.first);
		rhs.push_back(pair.second);
	}

	std::vector<IB_Column*> columns = read_columns(folder + "/columns", 0.5);
	std::cout << "columns read" << std::endl;


	return ISUD_Base(tasksNames, rhs, columns);

}

// Get problem from folder
ISUD_Base get_problem(std::string folder) {
	std::vector<std::pair<std::string, int>> tasks = read_tasks(folder + "/rhs.in");
	std::cout << "tasks read" << std::endl;
	std::vector<std::string> tasksNames;
	std::vector<int> rhs;
	for (auto pair : tasks) {
		tasksNames.push_back(pair.first);
		rhs.push_back(pair.second);
	}

	std::vector<IB_Column*> columns = read_columns(folder + "/columns");
	std::cout << "columns read" << std::endl;
	std::set<std::string> solution = read_initial_solution(folder + "/solution.out");
	std::cout << "solution read" << std::endl;
	std::map<std::string, int> contribs;
	std::set<std::string> done;
	int n = 0;
	for (int i = 0; i < columns.size(); i++) {
		if (solution.count(columns[i]->getName())) {
			done.insert(columns[i]->getName());
			n += 1;
			columns[i]->setInCurrentSolution();
			for (auto pair : columns[i]->getContribs()) {
				if (contribs.find(pair.first) != contribs.end()) {
					contribs[pair.first] += pair.second;
				}
				else {
					contribs[pair.first] = pair.second;
				}
			}
		}
		else {
			columns[i]->setOutCurrentSolution();
		}

		columns[i]->reorderContribs(tasksNames);
	}

	for (auto column : solution) {
		if (!done.count(column)) {
			std::cout << column << std::endl;
		}
	}


	std::cout << n << " " << solution.size() << std::endl;


	return ISUD_Base(tasksNames, rhs, columns);

}