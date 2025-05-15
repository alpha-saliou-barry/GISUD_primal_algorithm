#include "Problems.h"

// Generate all instances
void generateProblems() {
    std::vector<std::string> initial_folders = { "/home/barralph/Documents/GISUD_instances/prep/NW_320_0.34_0.33_0.33_1_day_16_day_23",
"/home/barralph/Documents/GISUD_instances/prep/NW_320_0.34_0.33_0.33_1_day_1_day_9",
"/home/barralph/Documents/GISUD_instances/prep/NW_757_0.34_0.33_0.33_1_day_12_day_22",
"/home/barralph/Documents/GISUD_instances/prep/NW_757_0.34_0.33_0.33_1_day_21_day_31" };

    std::vector<std::string> out_directories = { "/home/barralph/Documents/GISUD_instances_2/final/NW_320_0.34_0.33_0.33_1_day_16_day_23",
    "/home/barralph/Documents/GISUD_instances_2/final/NW_320_0.34_0.33_0.33_1_day_1_day_9",
    "/home/barralph/Documents/GISUD_instances_2/final/NW_757_0.34_0.33_0.33_1_day_12_day_22",
    "/home/barralph/Documents/GISUD_instances_2/final/NW_757_0.34_0.33_0.33_1_day_21_day_31" };

    //initial_folders = {"/home/barralph/Documents/GISUD_instances/prep/NW_757_0.34_0.33_0.33_1_day_21_day_31"};
    //out_directories = {"/home/barralph/Documents/GISUD_instances/final/NW_757_0.34_0.33_0.33_1_day_21_day_31"};
    std::vector<bool> changeRhs = { false, false, false, false };

    for (int k = 0; k < initial_folders.size(); k++) {
        std::string initial_folder = initial_folders[k];
        
        {
            ISUD_Base problem = generateProblemFromMps(initial_folder + "/columns.txt",
                initial_folder + "/rhs.txt", initial_folder + "/initial.txt", initial_folder + "/fixed_cost.txt");
            double max_cost = 0;
            for (int i = 0; i < problem.columns_.size(); i++) {
                if (problem.columns_[i]->getCost() > max_cost) {
                    max_cost = problem.columns_[i]->getCost();
                }
            }

            // Rajout des colonnes artificielles pour assurer la réalisabilité.
            std::set<int> acIndices;
            /*
            for (auto column : createInitialSolution(problem.tasks_, problem.rhs_, max_cost)) {
                problem.addColumn(column);
                acIndices.insert(problem.columns_.size() - 1);
            }*/

            std::cout << "Number of initial tasks : " << problem.tasks_.size() << std::endl;
            std::cout << "Number of columns : " << problem.columns_.size() << std::endl;
            std::map<int, int> rhs_numbers;
            for (int i = 0; i < problem.tasks_.size(); i++) {
                int rhs = problem.rhs_[i];
                if (rhs_numbers.find(rhs) == rhs_numbers.end()) {
                    rhs_numbers[rhs] = 0;
                }

                rhs_numbers[rhs] += 1;
            }

            for (auto pair : rhs_numbers) {
                std::cout << "Number of b_i = " << pair.first << " : " << pair.second << std::endl;
            }

            std::vector<int> initialSolution;
            for (int i = 0; i < problem.columns_.size(); i++) {
                initialSolution.push_back(problem.columns_[i]->isInCurrentSolution() ? 1 : 0);
            }

            ISUD isud(&problem, false, false, false);
            isud.solve(initial_folder);
            std::vector<int> solution = isud.getCurrentSolution();
            //std::cout << "Cplex objective : " << objective << std::endl;

            int n_ac_cols_in_solution = 0;
            int n_cols_in_solution = 0;
            for (int i = 0; i < solution.size(); i++) {
                if (solution[i] > 1e-4) {
                    problem.columns_[i]->setInCurrentSolution();
                    if (acIndices.count(i)) {
                        n_ac_cols_in_solution++;
                    }
                    n_cols_in_solution++;
                }
                else {
                    problem.columns_[i]->setOutCurrentSolution();
                }
            }

            std::ofstream characteristics;
            characteristics.open(out_directories[k] + "/characteristics.txt");
            characteristics << "Number of columns : " << problem.columns_.size() << std::endl;
            characteristics << "Number of tasks : " << problem.rhs_.size() << std::endl;
            characteristics << "Number of columns in solution : " << n_cols_in_solution << std::endl;
            characteristics << "Number of artificial variables in solution.: " << n_ac_cols_in_solution << std::endl;
            characteristics.close();

            computeInitialSolution(solution, 0.8, problem, out_directories[k] + "/perturbation_0.2");
            computeInitialSolution(solution, 0.7, problem, out_directories[k] + "/perturbation_0.3");
            computeInitialSolution(solution, 0.6, problem, out_directories[k] + "/perturbation_0.4");
        }
    }
}

// Export problem to a file in the GISUD format
// "problem" is the problem
// "columns_file" is the path to the columns file
// "rhs_file" is the path to the rhs file
// "initial_solution_file" is the path of the initial solution file
// "fixed_cost_path" is the path of the fixed cost file

void export_problem(ISUD_Base* problem, std::string columns_file, std::string rhs_file, std::string initial_solution_file, std::string fixed_cost_file) {
    std::ofstream file, rfile, isfile, fcfile;
    file.open(columns_file);
    rfile.open(rhs_file);

    if (fixed_cost_file != "") {
        fcfile.open(fixed_cost_file);
        fcfile << problem->fixed_cost_;
        fcfile.close();
    }

    if (initial_solution_file != "") {
        isfile.open(initial_solution_file);
    }

    std::map<std::string, int> tasksIndices;
    for (int i = 0; i < problem->tasks_.size(); i++) {
        tasksIndices[problem->tasks_[i]] = i;
    }

    file << problem->tasks_.size() << " " << problem->columns_.size() << std::endl;
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
    for (int i = 0; i < problem->rhs_.size(); i++) {
        rfile << problem->rhs_[i];
        if (i < problem->rhs_.size() - 1) {
            rfile << " ";
        }
    }

    rfile.close();

    if (initial_solution_file != "") {
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

}

// Construct ISUD_Base from paths
// "mps_file" is the mps_file path
// Compute initial solution for problem "problem"
// "solution" is the optimal solution of problem "problem"
// "perturbation_percent" is the percentage of remaining columns in new initial solution
// Out dir is the output directory where to put new problems with new initial solutions
int computeInitialSolution(std::vector<int> solution, double perturbation_percent, ISUD_Base& problem, std::string out_dir) {

    ISUD_Base initial_problem(problem);

    for (int i = 0; i < 4; i++) {
        perturbateFinalSolution(&problem, perturbation_percent);
        system(("mkdir -p " + out_dir + "/instance_" + std::to_string(i)).c_str());
        export_problem(&problem, out_dir + "/instance_" + std::to_string(i) + "/columns.txt",
            out_dir + "/instance_" + std::to_string(i) + "/rhs.txt",
            out_dir + "/instance_" + std::to_string(i) + "/initial.txt",
            out_dir + "/instance_" + std::to_string(i) + "/fixed_cost.txt");
        problem = initial_problem;
        for (int i = 0; i < solution.size(); i++) {
            if (solution[i] > 1e-4) {
                problem.columns_[i]->setInCurrentSolution();
            }
            else {
                problem.columns_[i]->setOutCurrentSolution();
            }
        }


    }

}

// "rhs_file" is the rhs_file path
// "initial_solution_file" is the initial solution path
// "fixed_cost_file" is the fixed cost file
ISUD_Base generateProblemFromMps(std::string mps_file, std::string rhs_file, std::string initial_solution_file, std::string fixed_cost_file) {
    std::ifstream file, rfile, isfile;
    if (rhs_file != "") {
        rfile.open(rhs_file);
    }
    std::cout << mps_file << std::endl;
    file.open(mps_file);
    std::string nRowsS, nColsS;
    file >> nRowsS;
    file >> nColsS;
    std::cout << nRowsS << " " << nColsS << std::endl;
    int nRows = std::stoi(nRowsS);
    std::vector<std::string> tasks;
    std::vector<int> rhs;
    for (int i = 0; i < nRows; i++) {
        tasks.push_back("task_" + std::to_string(i));
    }

    for (int i = 0; i < nRows; i++) {
        std::string rhs_;

        if (rhs_file != "") {
            rfile >> rhs_;
        }
        else {
            rhs_ = "1";
        }

        rhs.push_back(stoi(rhs_));
    }

    std::cout << rhs.size() << " " << tasks.size() << std::endl;
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
        if (initial_solution_file == "") {
            if (contribs.size() > 0 && ((rhs_vector - matrix(Eigen::seq(0, Eigen::last), columns.size() - 1)).array() >= 0).all()) {
                rhs_vector -= matrix(Eigen::seq(0, Eigen::last), columns.size() - 1);
                columns[columns.size() - 1]->setInCurrentSolution();
            }
            else {
                columns[columns.size() - 1]->setOutCurrentSolution();
            }
        }

    }

    file.close();

    int n_artificial_columns = 0;

    std::vector<int> remaining_rhs;
    for (int i = 0; i < rhs_vector.size(); i++) {
        remaining_rhs.push_back(rhs_vector(i));
    }

    std::cout << "Computation of artificial columns." << std::endl;

    if (initial_solution_file == "") {
        if (true) {
            for (auto column : addArtificialColumns(tasks, remaining_rhs, max_cost)) {
                columns.push_back(column);
                n_artificial_columns++;
            }
        }

    }
    else {
        isfile.open(initial_solution_file);
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
                std::cout << "Deletion of task : " << tasks[i] << std::endl;
                for (auto column : columns) {
                    column->deleteContrib(tasks[i]);
                }
            }
        }

        tasks = finalTasks;
        rhs = finalRhs;
    }


    //std::cout << rhs_vector.transpose() << std::endl;

    std::cout << n_artificial_columns << " artificial columns." << std::endl;

    std::vector<IB_Column*> final_columns;
    for (auto column : columns) {
        if (column->getContribs().size() > 0 || column->isInCurrentSolution()) {
            column->reorderContribs(tasks);
            final_columns.push_back(column);
        }
    }

    double fixed_cost = 0;
    if (fixed_cost_file != "") {
        std::ifstream fixed_cost_f;
        fixed_cost_f.open(fixed_cost_file);
        std::string fixed_cost_s;
        fixed_cost_f >> fixed_cost_s;
        fixed_cost = stod(fixed_cost_s);
        fixed_cost_f.close();
        std::cout << "Fixed cost : " << fixed_cost << std::endl;
    }

    return ISUD_Base(tasks, rhs, final_columns, fixed_cost);
}