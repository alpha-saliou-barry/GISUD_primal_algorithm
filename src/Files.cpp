#include "Files.h"

// Generate ISUD_Base problem from GISUD format
// Path of columns in "columns_path"
// Path of right-hand sides in "rhs_path"
// Path of initial solution in "initial_solution_path"
// Path of fixed cost file in "fixed_cost_path"

ISUD_Base constructISUDProblem(std::string columns_path, std::string rhs_path, std::string initial_solution_path, std::string fixed_cost_path) {
	

	std::ifstream file, rfile, isfile;
	if (rhs_path != "") {
		rfile.open(rhs_path);
	}
	//std::cout << mps_file << std::endl;

	file.open(columns_path);
	std::string nRowsS, nColsS;
	file >> nRowsS;
	file >> nColsS;
	// std::cout << nRowsS << " " << nColsS << std::endl;
	int nRows = std::stoi(nRowsS);
	std::vector<std::string> tasks;
	std::vector<int> rhs;
	for (int i = 0; i < nRows; i++) {
		tasks.push_back("task_" + std::to_string(i));
	}

	for (int i = 0; i < nRows; i++) {
		std::string rhs_;

		if (rhs_path != "") {
			rfile >> rhs_;
		}
		else {
			rhs_ = "1";
		}

		rhs.push_back(stoi(rhs_));
	}

	// std::cout << rhs.size() << " " << tasks.size() << std::endl;
	std::vector<IB_Column*> columns;
	Eigen::MatrixXi matrix(nRows, std::stoi(nColsS));
	matrix = Eigen::MatrixXi::Zero(nRows, std::stoi(nColsS));
	Eigen::VectorXi rhs_vector = Eigen::Map<Eigen::VectorXi>(rhs.data(), rhs.size());
	double max_cost = 0;
	while (file.good() && columns.size() < std::stoi(nColsS)) {

		std::string costS;
		std::string nNonZerosS;

		file >> costS;
		file >> nNonZerosS;
		std::string colName = "column_" + std::to_string(columns.size());

		double cost = std::stod(costS);
		std::vector<std::pair<std::string, int>> contribs;
		std::set<std::string> uniqueContribs;
		for (int i = 0; i < std::stoi(nNonZerosS); i++) {
			std::string rowContrib;
			file >> rowContrib;

			if (!uniqueContribs.count(tasks[((int)stoi(rowContrib))])) {
				contribs.push_back(std::pair<std::string, int>(tasks[((int)stoi(rowContrib))], 1));
				uniqueContribs.insert(tasks[((int)stoi(rowContrib))]);
				matrix(((int)stoi(rowContrib)), columns.size()) = 1;
			}
		}


		columns.push_back(new IB_Column(colName, contribs));
		if (cost > max_cost) {
			max_cost = cost;
		}
		columns[columns.size() - 1]->setCost(cost);
		columns[columns.size() - 1]->setOutCurrentSolution();
	}

	file.close();

	int n_artificial_columns = 0;

	std::vector<int> remaining_rhs;
	for (int i = 0; i < rhs_vector.size(); i++) {
		remaining_rhs.push_back(rhs_vector(i));
	}

	std::cout << "Recomputation of artificial columns." << std::endl;

	// Initial solution
	isfile.open(initial_solution_path);
	int n_pos_cols = 0;
	std::set<int> column_pos_nos;
	std::map<int, int> duplicates;

	while (isfile.good()) {
		std::string column_no;
		isfile >> column_no;
		columns[stoi(column_no)]->setInCurrentSolution();
		if (column_pos_nos.count(stoi(column_no))) {
			columns.push_back(new IB_Column(*columns[stoi(column_no)]));
			columns[columns.size() - 1]->setName(columns[columns.size() - 1]->getName() + "_DUPLICATE_" + std::to_string(duplicates[stoi(column_no)]));
			duplicates[stoi(column_no)] += 1;
			columns[columns.size() - 1]->setInCurrentSolution();
		}
		else {
			column_pos_nos.insert(stoi(column_no));
			duplicates[stoi(column_no)] = 1;
		}

		rhs_vector -= matrix(Eigen::seq(0, Eigen::last), stoi(column_no));
		n_pos_cols += 1;
	}

	std::cout << n_pos_cols << " positive columns." << std::endl;


	std::vector<std::string> finalTasks;
	std::vector<int> finalRhs;
	for (int i = 0; i < rhs_vector.size(); i++) {
		if (rhs_vector(i) == 0) {
			finalTasks.push_back(tasks[i]);
			finalRhs.push_back(rhs[i]);
		}
		else {
			std::cout << "Task deletion : " << tasks[i] << std::endl;
			for (auto column : columns) {
				column->deleteContrib(tasks[i]);
			}
		}
	}

	tasks = finalTasks;
	rhs = finalRhs;


	//std::cout << rhs_vector.transpose() << std::endl;

	//std::cout << n_artificial_columns << " artificial columns." << std::endl;

	std::vector<IB_Column*> final_columns;
	for (auto column : columns) {
		if (column->getContribs().size() > 0 || column->isInCurrentSolution()) {
			column->reorderContribs(tasks);
			final_columns.push_back(column);
		}
	}

	double fixed_cost = 0;
	if (fixed_cost_path != "") {
		std::ifstream fixed_cost_f;
		fixed_cost_f.open(fixed_cost_path);
		std::string fixed_cost_s;
		fixed_cost_f >> fixed_cost_s;
		fixed_cost = stod(fixed_cost_s);
		fixed_cost_f.close();
		std::cout << "Cout fixe : " << fixed_cost << std::endl;
	}

	return ISUD_Base(tasks, rhs, final_columns, fixed_cost);
}


// Export ISUD_Base problem to GISUD format
// Put columns in "columns_path" file
// Put right-hand-sides in "rhs_path" file
// Put fixed cost in "fixed_cost_path"
void export_problem(ISUD_Base* problem, std::string columns_path, std::string rhs_path, std::string initial_solution_path, std::string fixed_cost_path) {
	std::ofstream file, rfile, isfile, fcfile;
	file.open(columns_path);
	rfile.open(rhs_path);

	if (fixed_cost_path != "") {
		fcfile.open(fixed_cost_path);
		fcfile << problem->fixed_cost_;
		fcfile.close();
	}

	isfile.open(initial_solution_path);

	std::map<std::string, int> tasksIndices;
	for (int i = 0; i < problem->tasks_.size(); i++) {
		tasksIndices[problem->tasks_[i]] = i;
	}

	file << problem->tasks_.size() << " " << problem->columns_.size() << std::endl;

	// Columns file
	for (int i = 0; i < problem->columns_.size(); i++) {
		file << problem->columns_[i]->getCost() << " ";
		file << problem->columns_[i]->getContribs().size();

		for (std::string taskName : problem->columns_[i]->getContribsOrder()) {
			file << " " << tasksIndices[taskName];
		}

		if (i < problem->columns_.size() - 1) {
			file << std::endl;
		}
	}

	file.close();

	// Rhs file
	for (int i = 0; i < problem->rhs_.size(); i++) {
		rfile << problem->rhs_[i];
		if (i < problem->rhs_.size() - 1) {
			rfile << " ";
		}
	}

	rfile.close();

	// Initial solution file
	std::vector<int> initialSolution;
		for (int i = 0; i < problem->columns_.size(); i++) {
			if (problem->columns_[i]->isInCurrentSolution()) {
				initialSolution.push_back(i);
			}
		}

		for (int i = 0; i < initialSolution.size(); i++) {
			isfile << initialSolution[i];
			if (i < initialSolution.size() - 1) {
				isfile << " ";
			}
		}

		isfile.close();

}