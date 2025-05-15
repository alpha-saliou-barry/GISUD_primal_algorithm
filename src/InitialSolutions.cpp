#include "InitialSolutions.h"

// Pertubate final solution of problem "problem" at a given percetage "percentage".
// "percentage" is the percentage of columns that remains in final solution.
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
        // Choose randomly 2 columns
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
            std::cout << "Problem in the covered tasks" << std::endl;
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
        std::cout << "Number of columns : " << columns_in_solution.size() << std::endl;
    }

    for (auto column : columns_in_solution) {
        if (!problem->hasColumn(column->getName())) {
            problem->addColumn(column);
        }

        column->setInCurrentSolution();
    }
}

// Add artificial columns
// "tasks" is the problem tasks
// "remaining_rhs" is the remaining rhs to cover
// "max_cost" is the new cost given to each columns
std::vector<IB_Column*> addArtificialColumns(std::vector<std::string> tasks, std::vector<int> remaining_rhs, double max_cost) {
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

