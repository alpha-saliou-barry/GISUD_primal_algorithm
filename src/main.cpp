// ConsoleApplication1.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include "ISUD.h"
#include <string>
#include <ilcplex/ilocplex.h>
#include <cstring>
#include "CplexMIP.h"
#include <Eigen/Dense>
#include "IB_CompatibilityChecker.h"
#include <numeric>
#include "ProblemReader.h"



template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T>& v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    std::stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}

// Modify "percentage" of columns of the solution of problem
void perturbateFinalSolution(ISUD_Base* problem, double percentage = 0.5) {
    double max_cost = 0;
    std::vector<IB_Column*> columns_in_solution;
    std::set<std::string> initial_solution_columns;
    for (auto col : problem->columns_) {
        if (col->getCost() >= max_cost) {
            max_cost = col->getCost();
        }

        if (col->isInCurrentSolution()) {
            columns_in_solution.push_back(col);
            initial_solution_columns.insert(col->getName());
        }
    }

    std::cout << "Cout maximal : " << max_cost << std::endl;
    double remaining_percent = 1;
    int ac_number = 1;
    while (remaining_percent > percentage) {
        // Choose randomly 2 columns
        int first_column = rand() % columns_in_solution.size();

        IB_Column* first_column_vector = columns_in_solution[first_column];
        columns_in_solution.erase(columns_in_solution.begin() + first_column);

        std::set<std::string> coveredTasks;
        for (auto contrib : first_column_vector->getContribs()) {
            coveredTasks.insert(contrib.first);
        }

        std::vector<int> remainingColumns;
        for (int i = 0; i < columns_in_solution.size();i++) {
            IB_Column* column = columns_in_solution[i];
            bool take = true;
            for (auto contrib : column->getContribs()) {
                if (coveredTasks.count(contrib.first)) {
                    take = false;
                }
            }

            if (take) {
                remainingColumns.push_back(i);
            }
        }

        int second_column = rand() % (remainingColumns.size());
        IB_Column* second_column_vector = columns_in_solution[remainingColumns[second_column]];
        columns_in_solution.erase(columns_in_solution.begin() + remainingColumns[second_column]);
        for (auto contrib : second_column_vector->getContribs()) {
            coveredTasks.insert(contrib.first);
        }
        
        // Choose two random tasks in the column
        int first_task = rand() % first_column_vector->getContribs().size();
        int second_task = rand() % second_column_vector->getContribs().size();

        int counter = 0;
        std::string first_task_string, second_task_string;
        std::vector<std::pair<std::string, int>> newContribs1, newContribs2;
        std::set<std::string> coveredTasks1, coveredTasks2;
        bool save = true;

        for (auto taskName : first_column_vector->getContribsOrder()) {
            if (save) {
                coveredTasks1.insert(taskName);
                newContribs1.push_back(std::pair<std::string, int>(taskName, first_column_vector->findContribution(taskName)));
            }
            else {
                coveredTasks2.insert(taskName);
                newContribs2.push_back(std::pair<std::string, int>(taskName, first_column_vector->findContribution(taskName)));
            }
            
            if (counter == first_task) {
                first_task_string = taskName;
                save = false;
            } 

            counter++;
        }

        counter = 0;

        save = true;
        for (auto taskName : second_column_vector->getContribsOrder()) {
            if (save) {
                coveredTasks2.insert(taskName);
                newContribs2.push_back(std::pair<std::string, int>(taskName, second_column_vector->findContribution(taskName)));
            }
            else {
                coveredTasks1.insert(taskName);
                newContribs1.push_back(std::pair<std::string, int>(taskName, first_column_vector->findContribution(taskName)));
            }

            if (counter == second_task) {
                second_task_string = taskName;
                save = false;
            }

            counter++;
        }

        IB_Column* new_first_column = NULL, *new_second_column = NULL;
        for (auto column : problem->columns_) {
            if (!column->isInCurrentSolution() && new_first_column == NULL && column->isSame(coveredTasks1)) {
                new_first_column = column;
                new_first_column->setInCurrentSolution();
                columns_in_solution.push_back(new_first_column);
            } else if(!column->isInCurrentSolution() && (new_second_column == NULL && column->isSame(coveredTasks2))) {
                new_second_column = column;
                new_second_column->setInCurrentSolution();
                columns_in_solution.push_back(new_second_column);
            }
        }

        if (new_first_column == NULL && newContribs1.size() > 0) {
            new_first_column = new IB_Column("ac_column_" + std::to_string(ac_number), newContribs1);
            new_first_column->setCost(max_cost);
            new_first_column->setInCurrentSolution();
            ac_number++;
            columns_in_solution.push_back(new_first_column);
            new_first_column->reorderContribs(problem->tasks_);
        }
        if (new_second_column == NULL && newContribs2.size() > 0) {
            new_second_column = new IB_Column("ac_column_" + std::to_string(ac_number), newContribs2);
            new_second_column->setCost(max_cost);
            new_second_column->setInCurrentSolution();
            ac_number++;
            columns_in_solution.push_back(new_second_column);
            new_second_column->reorderContribs(problem->tasks_);
        }
        std::set<std::string> dest;
        std::set_union(coveredTasks1.begin(), coveredTasks1.end(), coveredTasks2.begin(), coveredTasks2.end(), std::inserter(dest, dest.begin()));
        if (dest != coveredTasks) {
            std::cout << "Il y a un problème dans les taches couvertes." << std::endl;
            throw std::exception();
        }

        first_column_vector->setOutCurrentSolution();
        second_column_vector->setOutCurrentSolution();
        int n_good_columns = 0;
        for (auto col : columns_in_solution) {
            if (initial_solution_columns.count(col->getName())) {
                n_good_columns += 1;
            }
        }

        remaining_percent = ((float) n_good_columns) / initial_solution_columns.size();
        std::cout << "Nouveau pourcentage : " << remaining_percent << std::endl;
        std::cout << "Nombre de colonnes : " << columns_in_solution.size() << std::endl;
    }

    for (auto column : columns_in_solution) {
        if (!problem->hasColumn(column->getName())) {
            problem->addColumn(column);
        }

        column->setInCurrentSolution();
    }
}

// Compute initial solutions
std::vector<IB_Column*> createInitialSolution(std::vector<std::string> tasks, std::vector<int> remaining_rhs, double max_cost) {
    std::vector<IB_Column*> result;
    int n_tasks = 0;
    int current_task = -1;
    int nb_covered_rhs = 0;
    for (int i = 0; i < remaining_rhs.size();i++) {
        int rhs = remaining_rhs[i];
        if (rhs > 0 && current_task == -1) {
            current_task = i;
        }

        n_tasks += rhs;
    }

    for (int i = 0; i < n_tasks; i++) {
        std::vector<std::pair<std::string, int>> contribs;
        contribs.push_back(std::pair<std::string, int>(tasks[current_task], 1));
        result.push_back(new IB_Column("cartificial_" + std::to_string(i), contribs));
        result[result.size() - 1]->setInCurrentSolution();
        result[result.size() - 1]->setCost(max_cost);
        result[result.size() - 1]->reorderContribs(tasks);
        nb_covered_rhs += 1;
        if (nb_covered_rhs >= remaining_rhs[current_task]) {
            nb_covered_rhs = 0;
            int j = current_task + 1;
            while (j < remaining_rhs.size() && remaining_rhs[j] == 0) {
                j++;
            }

            if (j < remaining_rhs.size()) {
                current_task = j;
            }
            else {
                break;
            }
        }
    }

    return result;
}

// Export problem to a file in the GISUD format
void export_problem(ISUD_Base* problem, std::string columns_file, std::string rhs_file, std::string initial_solution_file, std::string fixed_cost_file = "") {
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
    for (int i = 0; i < problem->columns_.size();i++) {
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

// Read a problem from the GISUD format
ISUD_Base generateProblemFromMps(std::string mps_file, std::string rhs_file, std::string initial_solution_file = "", std::string fixed_cost_file = "") {
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

            if(!uniqueContribs.count(tasks[((int) stoi(rowContrib))])) {
                contribs.push_back(std::pair<std::string, int>(tasks[((int) stoi(rowContrib))], 1));
                uniqueContribs.insert(tasks[((int) stoi(rowContrib))]);
                matrix(((int) stoi(rowContrib)), columns.size()) = 1;
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

    std::cout << "Calcul des colonnes artificiels." << std::endl;

    if (initial_solution_file == "") {
        if (true) {
            for (auto column : createInitialSolution(tasks, remaining_rhs, max_cost)) {
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
        while(isfile.good()) {
            std::string column_no;
            isfile >> column_no;
            columns[stoi(column_no)]->setInCurrentSolution();
            if(column_pos_nos.count(stoi(column_no))) {
                columns.push_back(new IB_Column(*columns[stoi(column_no)]));
                columns[columns.size()-1]->setName(columns[columns.size()-1]->getName() + "_DUPLICATE_" + std::to_string(duplicates[stoi(column_no)]));
                duplicates[stoi(column_no)] += 1;
                columns[columns.size()-1]->setInCurrentSolution();
            } else {
                column_pos_nos.insert(stoi(column_no));
                duplicates[stoi(column_no)] = 1;
            }
            
            rhs_vector -= matrix(Eigen::seq(0, Eigen::last), stoi(column_no));
            n_pos_cols += 1;
        }

        std::cout << n_pos_cols << " colonnes positives." << std::endl;

        
        std::vector<std::string> finalTasks;
        std::vector<int> finalRhs;
        for(int i = 0; i < rhs_vector.size();i++) {
            if(rhs_vector(i) == 0) {
                finalTasks.push_back(tasks[i]);
                finalRhs.push_back(rhs[i]);
            } else {
                std::cout << "Suppression de la tache : " << tasks[i] << std::endl;
                for(auto column: columns) {
                    column->deleteContrib(tasks[i]);
                }
            }
        }

        tasks = finalTasks;
        rhs = finalRhs;
    }
    

    std::cout << rhs_vector.transpose() << std::endl;

    std::cout << n_artificial_columns << " colonnes artificielles." << std::endl;

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
        std::cout << "Cout fixe : " << fixed_cost << std::endl;
    }

    return ISUD_Base(tasks, rhs, final_columns, fixed_cost);
}


// Create problems with different initial solutions from the solution solution and the problem problem
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

// Change the RHS
void perturbateRhs(ISUD_Base* problem, float perturbation_percent = 0.25) {
    int first_rhs = 0;
    if (perturbation_percent < 1 - 1e-4) {
        first_rhs = rand() % ((int)(problem->rhs_.size() * (1 - perturbation_percent)));
    }


    int i = first_rhs;

    while (i < problem->rhs_.size() && i - first_rhs <= perturbation_percent * problem->rhs_.size()) {
        problem->rhs_[i] = 2;
        i += 1;
    }
}

// Convert problem in the ORLIB format to the GISUD format
void convertProblems() {
    std::vector < std::pair <std::string, std::string> > problems = { std::pair<std::string, std::string>("/home/barralph/Documents/GISUD_instances/initial/sppaa05/columns.txt", ""),
    std::pair<std::string, std::string>("/home/barralph/Documents/GISUD_instances/initial/NW_320/columns.txt", "/home/barralph/Documents/GISUD_instances/initial/NW_320/rhs.txt"),
    std::pair<std::string, std::string>("/home/barralph/Documents/GISUD_instances/initial/NW_757/columns.txt", "/home/barralph/Documents/GISUD_instances/initial/NW_757/rhs.txt") };

    std::vector<std::string> out_directories = { "/home/barralph/Documents/GISUD_instances/initial/sppaa05","/home/barralph/Documents/GISUD_instances/initial/NW_320",
        "/home/barralph/Documents/GISUD_instances/initial/NW_757" };
    for (int i = 0; i < problems.size();i++) {
        std::pair<std::string, std::string> pair = problems[i];
        ISUD_Base problem = generateProblemFromMps(pair.first,
            pair.second);
        export_problem(&problem, out_directories[i] + "/columns2.txt", out_directories[i] + "/rhs2.txt", "");
    }
}

// Compute the statistics of our problems
void getStatistics() {
    /*std::vector<std::string> initial_folders = { "/home/barralph/Documents/GISUD_instances/initial/sppaa01",
        "/home/barralph/Documents/GISUD_instances/initial/instance_757",
        "/home/barralph/Documents/GISUD_instances/initial/instance_320",
        "/home/barralph/Documents/GISUD_instances/prep/NW_320" };

    std::vector<std::string> problem_names = { "sppaa01", "instance_757", "instance_320", "NW_320" };

    std::ofstream out_file;
    out_file.open("stats.csv");
    out_file << "Instance;N columns;N tasks;Average non zeros" << std::endl;
    for (int k = 0; k < initial_folders.size(); k++) {
        std::string initial_folder = initial_folders[k];
        ISUD_Base problem = generateProblemFromMps(initial_folder + "/columns.txt",
            initial_folder + "/rhs.txt", "", initial_folder + "/fixed_cost.txt");

        int n_non_zeros = 0;
        for (auto column : problem.columns_) {
            n_non_zeros += column->getContribs().size();
        }

        out_file << problem_names[k] << ";" << problem.columns_.size() << ";" << problem.tasks_.size() << ";" << ((double)n_non_zeros) / problem.columns_.size() << std::endl;
    }

    out_file.close();*/

    std::vector<std::string> initial_folders = { "/home/barralph/Documents/GISUD_instances/final/NW_320_0.25_0.25_0.25_0.25_1_day_16_day_23/perturbation_0.2/instance_0",
    "/home/barralph/Documents/GISUD_instances/final/NW_320_0.25_0.25_0.25_0.25_1_day_1_day_9/perturbation_0.2/instance_0",
    "/home/barralph/Documents/GISUD_instances/final/NW_320_0.34_0.33_0.33_1_day_16_day_23/perturbation_0.2/instance_0",
    "/home/barralph/Documents/GISUD_instances/final/NW_320_0.34_0.33_0.33_1_day_1_day_9/perturbation_0.2/instance_0",
     "/home/barralph/Documents/GISUD_instances/final/NW_757_0.25_0.25_0.25_0.25_1_day_12_day_22/perturbation_0.2/instance_0",
     "/home/barralph/Documents/GISUD_instances/final/NW_757_0.25_0.25_0.25_0.25_1_day_1_day_11/perturbation_0.2/instance_0",
     "/home/barralph/Documents/GISUD_instances/final/NW_757_0.34_0.33_0.33_1_day_12_day_22/perturbation_0.2/instance_0",
     "/home/barralph/Documents/GISUD_instances/final/NW_757_0.34_0.33_0.33_1_day_1_day_11/perturbation_0.2/instance_0" };;

    std::vector<std::string> problem_names = { "320_3", "320_4", "320_1", "320_2", "757_3", "757_4", "757_1", "757_2"};
    std::ofstream out_file;
    out_file.open("rhs.txt");

    for (int k = 0; k < initial_folders.size(); k++) {
        std::string initial_folder = initial_folders[k];
        ISUD_Base problem = generateProblemFromMps(initial_folder + "/columns.txt",
            initial_folder + "/rhs.txt", "");

        int n_non_zeros = 0;
        for (auto column : problem.columns_) {
            n_non_zeros += column->getContribs().size();
        }

        out_file << std::endl << problem_names[k] << std::endl;
        out_file << "Nombre de lignes initial : " << problem.tasks_.size() << std::endl;
        //std::cout << "Nombre de lignes agrégées : " << problem_aggregated.tasks_.size() << std::endl;
        out_file << "Nombre de colonnes : " << problem.columns_.size() << std::endl;
        out_file << "Densite : " <<  ((double)n_non_zeros) / problem.columns_.size() << std::endl;

        std::map<int, int> rhs_numbers;
        for (int i = 0; i < problem.tasks_.size(); i++) {
            int rhs = problem.rhs_[i];
            if (rhs_numbers.find(rhs) == rhs_numbers.end()) {
                rhs_numbers[rhs] = 0;
            }

            rhs_numbers[rhs] += 1;
        }

        for (auto pair : rhs_numbers) {
            out_file << "Nombre de b_i = " << pair.first << " : " << pair.second << std::endl;
        }
    }

    out_file.close();

    std::cout << "Statistiques générées" << std::endl;
}

// Generate all instances
void generateProblems() {
        std::vector<std::string> initial_folders = { "/home/barralph/Documents/GISUD_instances/prep/NW_320_0.34_0.33_0.33_1_day_16_day_23",
    "/home/barralph/Documents/GISUD_instances/prep/NW_320_0.34_0.33_0.33_1_day_1_day_9",
    "/home/barralph/Documents/GISUD_instances/prep/NW_757_0.34_0.33_0.33_1_day_12_day_22",
    "/home/barralph/Documents/GISUD_instances/prep/NW_757_0.34_0.33_0.33_1_day_21_day_31"};

    std::vector<std::string> out_directories = { "/home/barralph/Documents/GISUD_instances_2/final/NW_320_0.34_0.33_0.33_1_day_16_day_23",
    "/home/barralph/Documents/GISUD_instances_2/final/NW_320_0.34_0.33_0.33_1_day_1_day_9",
    "/home/barralph/Documents/GISUD_instances_2/final/NW_757_0.34_0.33_0.33_1_day_12_day_22",
    "/home/barralph/Documents/GISUD_instances_2/final/NW_757_0.34_0.33_0.33_1_day_21_day_31"};
    
    //initial_folders = {"/home/barralph/Documents/GISUD_instances/prep/NW_757_0.34_0.33_0.33_1_day_21_day_31"};
    //out_directories = {"/home/barralph/Documents/GISUD_instances/final/NW_757_0.34_0.33_0.33_1_day_21_day_31"};
    std::vector<bool> changeRhs = { false, false, false, false };

    for (int k = 0; k < initial_folders.size();k++) {
        std::string initial_folder = initial_folders[k];
        if (changeRhs[k])  {
            int rhs_no = 1;
            for (auto perturbation : { 0.25, 0.5 }) {
                ISUD_Base problem = generateProblemFromMps(initial_folder + "/columns.txt",
                    initial_folder + "/rhs.txt", "", initial_folder + "/fixed_cost.txt");
                perturbateRhs(&problem, perturbation);
                double max_cost = 0;
                for (int i = 0; i < problem.columns_.size(); i++) {
                    if (problem.columns_[i]->getCost() > max_cost) {
                        max_cost = problem.columns_[i]->getCost();
                    }
                }

                std::set<int> acIndices;
                // Rajout des colonnes artificielles pour assurer la réalisabilité.
                for (auto column : createInitialSolution(problem.tasks_, problem.rhs_, max_cost)) {
                    problem.addColumn(column);
                    acIndices.insert(problem.columns_.size() - 1);
                }
                
                std::cout << "Nombre de lignes initial : " << problem.tasks_.size() << std::endl;
                std::cout << "Nombre de colonnes : " << problem.columns_.size() << std::endl;
                std::map<int, int> rhs_numbers;
                for (int i = 0; i < problem.tasks_.size(); i++) {
                    int rhs = problem.rhs_[i];
                    if (rhs_numbers.find(rhs) == rhs_numbers.end()) {
                        rhs_numbers[rhs] = 0;
                    }

                    rhs_numbers[rhs] += 1;
                }

                for (auto pair : rhs_numbers) {
                    std::cout << "Nombre de b_i = " << pair.first << " : " << pair.second << std::endl;
                }

                std::vector<int> initialSolution;
                CplexMIP cplex(&problem);
                std::vector<int> solution;
                double objective = cplex.solve(initialSolution, &solution);
                std::cout << "Cplex objective : " << objective << std::endl;


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
                characteristics.open("/home/barralph/Documents/GISUD_instances/characteristics_" + std::to_string(k) + "_" + std::to_string(perturbation) +".txt");
                characteristics << "Nombre de colonnes : " << problem.columns_.size() << std::endl;
                characteristics << "Nombre de taches : " << problem.rhs_.size() << std::endl;
                characteristics << "Nombre de colonnes dans la solution : " << n_cols_in_solution << std::endl;
                characteristics << "Nombre de colonnes de variables artificielles dans la solution : " << n_ac_cols_in_solution << std::endl;
                characteristics.close();

                computeInitialSolution(solution, 0.95, problem, out_directories[k] + "/rhs_"+ std::to_string(rhs_no) + "/perturbation_0.05");
                computeInitialSolution(solution, 0.975, problem, out_directories[k] + "/rhs_" + std::to_string(rhs_no) + "/perturbation_0.025");
                //computeInitialSolution(solution, 0.6, problem, out_directories[k] + "/rhs_" + std::to_string(rhs_no) + "/perturbation_0.4");

                rhs_no += 1;
            }
        }
        else {
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

            std::cout << "Nombre de lignes initial : " << problem.tasks_.size() << std::endl;
                std::cout << "Nombre de colonnes : " << problem.columns_.size() << std::endl;
            std::map<int, int> rhs_numbers;
            for (int i = 0; i < problem.tasks_.size(); i++) {
                int rhs = problem.rhs_[i];
                if (rhs_numbers.find(rhs) == rhs_numbers.end()) {
                    rhs_numbers[rhs] = 0;
                }

                rhs_numbers[rhs] += 1;
            }

            for (auto pair : rhs_numbers) {
                std::cout << "Nombre de b_i = " << pair.first << " : " << pair.second << std::endl;
            }

            std::vector<int> initialSolution;
            for(int i = 0; i < problem.columns_.size();i++) {
                initialSolution.push_back(problem.columns_[i]->isInCurrentSolution() ? 1:0);
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
            characteristics << "Nombre de colonnes : " << problem.columns_.size() << std::endl;
            characteristics << "Nombre de taches : " << problem.rhs_.size() << std::endl;
            characteristics << "Nombre de colonnes dans la solution : " << n_cols_in_solution << std::endl;
            characteristics << "Nombre de colonnes de variables artificielles dans la solution : " << n_ac_cols_in_solution << std::endl;
            characteristics.close();

            computeInitialSolution(solution, 0.8, problem, out_directories[k] + "/perturbation_0.2");
            computeInitialSolution(solution, 0.7, problem, out_directories[k] +  "/perturbation_0.3");
            computeInitialSolution(solution, 0.6, problem, out_directories[k] +  "/perturbation_0.4");
        }
    }
}

// Compute statistics of our problems
void computeCharacteristics() {
    std::vector<std::string> initialPaths = {
                                             "/home/barralph/Documents/GISUD_instances/final/instance_320/rhs_1",
                                             "/home/barralph/Documents/GISUD_instances/final/instance_320/rhs_2",
                                             "/home/barralph/Documents/GISUD_instances/final/instance_757/rhs_1",
                                             "/home/barralph/Documents/GISUD_instances/final/instance_757/rhs_2",
    "/home/barralph/Documents/GISUD_instances/final/NW_757",
    "/home/barralph/Documents/GISUD_instances/final/NW_320" };

    std::vector<std::string> paths = {  };
    for (int i = 0; i < initialPaths.size(); i++) {
        std::string path = initialPaths[i];
        for(auto perturbation : {"0.2"  }) {
            for (int j = 0; j < 1; j++) {
                paths.push_back(path + "/perturbation_" + perturbation + "/instance_" + std::to_string(j));
            }
        }
    }

    for (int i = 0; i < paths.size(); i++) {
        std::string path = paths[i];

        ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt", path + "/fixed_cost.txt");
        std::vector<int> solution;
        std::set<int> acIndices;
        double maxCost = 0;
        int n_non_zeros = 0;
        for (int j = 0; j < problem.columns_.size(); j++) {
            if (problem.columns_[j]->getCost() >= maxCost) {
                maxCost = problem.columns_[j]->getCost();
            }

            n_non_zeros += problem.columns_[j]->getContribs().size();
        }

        for (int j = 0; j < problem.columns_.size(); j++) {
            solution.push_back(problem.columns_[j]->isInCurrentSolution() ? 1 : 0 );
            if (problem.columns_[j]->isInCurrentSolution() && problem.columns_[j]->getCost() >= maxCost - 1e3) {
                acIndices.insert(j);
            }
        }

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
        std::map<int, int> rhs_numbers;
        for (int i = 0; i < problem.tasks_.size(); i++) {
            int rhs = problem.rhs_[i];
            if (rhs_numbers.find(rhs) == rhs_numbers.end()) {
                rhs_numbers[rhs] = 0;
            }

            rhs_numbers[rhs] += 1;
        }

        
        std::ofstream characteristics;
        characteristics.open(path + "/characteristics.txt");
        characteristics << "Nombre de colonnes : " << problem.columns_.size() << std::endl;
        characteristics << "Nombre de taches : " << problem.rhs_.size() << std::endl;
        characteristics << "Nombre de colonnes dans la solution : " << n_cols_in_solution << std::endl;
        characteristics << "Nombre de colonnes de variables artificielles dans la solution : " << n_ac_cols_in_solution << std::endl;
        characteristics << "Nombre moyen de non zeros par colonne : " << ((double)n_non_zeros) / problem.columns_.size() << std::endl;
        for (auto pair : rhs_numbers) {
          characteristics <<           "Nombre de b_i = " << pair.first << " : " << pair.second << std::endl;
        }
        characteristics.close();
    }
    
}

// Solve with CPLEX the problem
void calculateInstanceCplex() {
    std::vector<std::string> paths = { "/home/barralph/Documents/newInstances/instance_320/rhs_1531_936_305_158" };
    for (auto path : paths) {
        ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt");
        double max_cost = 0;
        for (auto column : problem.columns_) {
            if (column->getCost() >= max_cost) {
                max_cost = column->getCost();
            }
        }

        for (auto column : createInitialSolution(problem.tasks_, problem.rhs_, max_cost)) {
            problem.addColumn(column);
        }

        // On calcule le pre
        CplexMIP cplex(&problem);
        std::vector<int> currentSolution;
        std::vector<int> solution;
        cplex.solve(currentSolution, &solution, false, "", path + "/problem.mps");
        //std::rename((path + "/problem.pre").c_str(), (path + "/problem.mps").c_str());
        //cplex.solveFromFile(path + "/problem.mps", path + "/solution.txt");
    }
}


void calculatePreSolution() {
    std::vector<std::string> paths = { "/home/barralph/instances_gspp/NW_320" };
    for (auto path : paths) {
        CplexMIP cplex(NULL);
        cplex.solveFromFile(path + "/problem.mps", path + "/solution.txt");
    }
}

// Compute CPLEX solutions of our instances
int computeCplexSolutions() {
    std::vector<std::string> out_directories = {
        "/home/barralph/Documents/GISUD_instances/final/0.25_0.25_0.25_0.25_1_day_16_day_23",
    "/home/barralph/Documents/GISUD_instances/final/0.34_0.33_0.33_1_day_16_day_23"};
    //out_directories = {"/home/barralph/Documents/GISUD_instances/final/alphainstance"};
    std::vector<std::string> rhs_directories = { "" };
    std::vector<std::string> perturbations_levels = { "perturbation_0.2","perturbation_0.3", "perturbation_0.4"};
    std::vector<std::string> instances = { "instance_0","instance_1", "instance_2", "instance_3", "instance_4", "instance_5",
    "instance_6", "instance_7" };
    

    for (std::string out_directory : out_directories) {
        for (std::string rhs_directory : rhs_directories) {
            for (std::string perturbation_level : perturbations_levels) {
                for (std::string instance : instances) {
                    std::string path = out_directory + "/" + rhs_directory + "/"+ perturbation_level + "/"+ instance;
                    {
                    ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
                    ISUD isud(&problem, false, false, true);
                    std::vector<int> initialSolution = isud.getCurrentSolution();
                    //CplexMIP cplex(&problem);
                    //std::vector<int> solution;

                    //double objective = cplex.solve(initialSolution, &solution, false, path);
                    //std::ifstream outCplex;
                    //outCplex.open(path + "/cplex_solution.txt");
                    //std::string objective2, bound;
                    //outCplex >> objective2;
                    //outCplex >> bound;
                    //outCplex.close();
                    //double bestBound = std::stod(bound);


                    isud.solve(path);
                    problem.destroyColumns();
                    }
                    {
                    ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
                    ISUD isud(&problem, false, false, false);
                    std::vector<int> initialSolution = isud.getCurrentSolution();
                    //CplexMIP cplex(&problem);
                    //std::vector<int> solution;

                    //double objective = cplex.solve(initialSolution, &solution, false, path);
                    //std::ifstream outCplex;
                    //outCplex.open(path + "/cplex_solution.txt");
                    //std::string objective2, bound;
                    //outCplex >> objective2;
                    //outCplex >> bound;
                    //outCplex.close();
                    //double bestBound = std::stod(bound);


                    isud.solve(path);
                    problem.destroyColumns();
                    }
                    {
                        ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
                        ISUD isud(&problem, true, false, false);
                        std::vector<int> initialSolution = isud.getCurrentSolution();
                        //CplexMIP cplex(&problem);
                        //std::vector<int> solution;

                        //double objective = cplex.solve(initialSolution, &solution, false, path);
                        //std::ifstream outCplex;
                        //outCplex.open(path + "/cplex_solution.txt");
                        //std::string objective2, bound;
                        //outCplex >> objective2;
                        //outCplex >> bound;
                        //outCplex.close();
                        //double bestBound = std::stod(bound);


                        isud.solve(path);
                        problem.destroyColumns();
                    }
                }
            }
        }
    }
}

// Generate the problem with its initial solution
int computeInitialSolution2(std::string path_) {
    std::vector<std::string> paths = { path_ };
    for (auto path : paths) {
        ISUD_Base problem = get_problem(path);
        CplexMIP cplex(&problem);
        ISUD isud(&problem, true, false);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        std::vector<int> solution;
        double objective = cplex.solve(initialSolution, &solution, false);
        std::cout << "Objectif : " << objective << std::endl;
        for (int i = 0; i < solution.size(); i++) {
            if (solution[i] == 1) {
                problem.columns_[i]->setInCurrentSolution();
            }
            else {
                problem.columns_[i]->setOutCurrentSolution();
            }
        }

        perturbateFinalSolution(&problem, 0.7);
        export_problem(&problem, path + "/columns_alpha.txt", path + "/rhs_alpha.txt", path + "/initial_solution_alpha.txt");
        problem.destroyColumns();
    }
}


// Compute CPLEX solutions of our instances
int computeCplexSolutions2(std::string path_) {
    std::vector<std::string> paths = { path_ };
    for (auto path : paths) {
        ISUD_Base problem = generateProblemFromMps(path + "/columns_alpha.txt", path + "/rhs_alpha.txt", path + "/initial_solution_alpha.txt");
        ISUD isud(&problem, true, false);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        CplexMIP cplex(&problem);
        std::vector<int> solution;

        double objective = cplex.solve(initialSolution, &solution, false, path);
        std::cout << "Objectif : " << objective << std::endl;
        problem.destroyColumns();
    }
}


// Compute GISUD solutions of our instances with disaggregation
int computeISUDSolutions2(std::string path_) {
    std::vector<std::string> paths = { path_ };
    for (auto path : paths) {
        ISUD_Base problem = generateProblemFromMps(path + "/columns_alpha.txt", path + "/rhs_alpha.txt", path + "/initial_solution_alpha.txt");
        ISUD isud(&problem, true, false);
        std::ifstream outCplex;
        outCplex.open(path + "/cplex_solution.txt");
        std::string objective, bound;
        outCplex >> objective;
        outCplex >> bound;
        outCplex.close();
        double bestBound = std::stod(bound);
        std::cout << bestBound << std::endl;

        isud.solve(path);
        problem.destroyColumns();
    }
}

// Compute GISUD solutions without disaggregation and column addition
int computeISUDSolutions() {
    std::vector<std::string> out_directories = { "/home/barralph/Documents/GISUD_instances/prep/ninstance" };

    std::vector<std::string> rhs_directories = { ""};
    std::vector<std::string> perturbations_levels = { "perturbation_0.2", "perturbation_0.3", "perturbation_0.4"};
    perturbations_levels = {""};
    std::vector<std::string> instances = { "instance_0", "instance_1", "instance_2", "instance_3", "instance_4", "instance_5", "instance_6", "instance_7" };
     instances = {""};

    for (std::string out_directory : out_directories) {
        for (std::string rhs_directory : rhs_directories) {
            for (std::string perturbation_level : perturbations_levels) {
                for (std::string instance : instances) {
                    std::string path = out_directory + "/" + rhs_directory + "/" + perturbation_level + "/" + instance;
                    ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt", path + "/fixed_cost.txt");
                    ISUD isud(&problem, false, false);
                    double bestBound = 0;

                    std::vector<int> initialSolution = isud.getCurrentSolution();
                    std::cout << bestBound << std::endl;
                    isud.solve(path);
                    problem.destroyColumns();
                }
            }
        }
    }
}


// Main function
int main(int argc, char* argv[])
{
  //modifyRhs();
  //return 0;
  //modifyRhs("/home/barralph/Documents/TestsArticles/vcsMB50_1", "/home/barralph/Documents/TestsArticles/vcsMB50_2");
  //generateProblems();
    //computeCplexSolutions();
  //return 0;
    std::string path = argv[1];
    char* cplex = argv[2];
    if(strcmp(cplex, "0") == 0) {
        char* rc = argv[3];
        char* dis = argv[4];
        char* compete = argv[5];
        ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
        ISUD isud(&problem, strcmp(rc, "1") == 0, false, strcmp(dis, "1") == 0, strcmp(compete, "1") == 0);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        isud.solve(path);
        problem.destroyColumns();
    } else {
      ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt","",
                        path + "/fixed_cost.txt");
        ISUD isud(&problem, false, false, true);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        CplexMIP cplex(&problem);
        std::vector<int> solution;

        double objective = cplex.solve(initialSolution, &solution, false, path);
        problem.destroyColumns();
    }

    return 0;
}
