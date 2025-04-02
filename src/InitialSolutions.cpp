#include "InitialSolutions.h"


// Perturbate the optimal solution of the problem problem at a given percentage
void perturbateFinalSolution(ISUD_Base* problem, double percentage) {
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

    std::cout << "Maximal cost : " << max_cost << std::endl;
    double remaining_percent = 1;
    int ac_number = 1;
    while (remaining_percent > percentage) {
        // Two random columns are chosen.
        int first_column = rand() % columns_in_solution.size();

        IB_Column* first_column_vector = columns_in_solution[first_column];
        columns_in_solution.erase(columns_in_solution.begin() + first_column);

        std::set<std::string> coveredTasks;
        for (auto contrib : first_column_vector->getContribs()) {
            coveredTasks.insert(contrib.first);
        }

        std::vector<int> remainingColumns;
        for (int i = 0; i < columns_in_solution.size(); i++) {
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

        // Two random tasks are chosen. We connect these two tasks.
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

        IB_Column* new_first_column = NULL, * new_second_column = NULL;
        for (auto column : problem->columns_) {
            if (!column->isInCurrentSolution() && new_first_column == NULL && column->isSame(coveredTasks1)) {
                new_first_column = column;
                new_first_column->setInCurrentSolution();
                columns_in_solution.push_back(new_first_column);
            }
            else if (!column->isInCurrentSolution() && (new_second_column == NULL && column->isSame(coveredTasks2))) {
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
            std::cout << "Problem in covered tasks." << std::endl;
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

        remaining_percent = ((float)n_good_columns) / initial_solution_columns.size();
        std::cout << "New percentage : " << remaining_percent << std::endl;
        std::cout << "Columns number : " << columns_in_solution.size() << std::endl;
    }

    for (auto column : columns_in_solution) {
        if (!problem->hasColumn(column->getName())) {
            problem->addColumn(column);
        }

        column->setInCurrentSolution();
    }
}

// Return artificial columns of cost "max_cost" to add to a problem on tasks "tasks" with a remaining rhs "remaining_rhs"
std::vector<IB_Column*> getArtificialColumns(std::vector<std::string> tasks, std::vector<int> remaining_rhs, double max_cost) {
    std::vector<IB_Column*> result;
    int n_tasks = 0;
    int current_task = -1;
    int nb_covered_rhs = 0;
    for (int i = 0; i < remaining_rhs.size(); i++) {
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


// Compute instances by perturbation from a ISUD problem with optimal solution solution 
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

// Compute "n_instances" instances by perturbation from a ISUD problem with optimal solution "solution" to output directory "out_dir"
int computeInitialSolution(std::vector<int> solution, double perturbation_percent, ISUD_Base& problem, std::string out_dir, int n_instances) {
    ISUD_Base initial_problem(problem);

    for (int i = 0; i < n_instances; i++) {
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

// Generate problems from initial_folders to out_directories by perturbing the final solution
void generateProblems(std::vector<std::string> initial_folders, std::vector<std::string> out_directories) {
    for (int k = 0; k < initial_folders.size(); k++) {
        std::string initial_folder = initial_folders[k];
        
        {
            ISUD_Base problem = constructISUDProblem(initial_folder + "/columns.txt",
                initial_folder + "/rhs.txt", initial_folder + "/initial.txt", initial_folder + "/fixed_cost.txt");
            double max_cost = 0;
            for (int i = 0; i < problem.columns_.size(); i++) {
                if (problem.columns_[i]->getCost() > max_cost) {
                    max_cost = problem.columns_[i]->getCost();
                }
            }

            // Add artificial columns
            std::set<int> acIndices;

            std::cout << "Initial rows number : " << problem.tasks_.size() << std::endl;
            std::cout << "Columns number : " << problem.columns_.size() << std::endl;
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

            computeInitialSolution(solution, 0.8, problem, out_directories[k] + "/perturbation_0.2", 4);
            computeInitialSolution(solution, 0.7, problem, out_directories[k] + "/perturbation_0.3", 4);
            computeInitialSolution(solution, 0.6, problem, out_directories[k] + "/perturbation_0.4", 4);
        }
    }
}