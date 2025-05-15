#include "Statistics.h"

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

    std::vector<std::string> problem_names = { "320_3", "320_4", "320_1", "320_2", "757_3", "757_4", "757_1", "757_2" };
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
        out_file << "Number of initial tasks : " << problem.tasks_.size() << std::endl;
        //std::cout << "Nombre de lignes agrégées : " << problem_aggregated.tasks_.size() << std::endl;
        out_file << "Columns number : " << problem.columns_.size() << std::endl;
        out_file << "Density : " << ((double)n_non_zeros) / problem.columns_.size() << std::endl;

        std::map<int, int> rhs_numbers;
        for (int i = 0; i < problem.tasks_.size(); i++) {
            int rhs = problem.rhs_[i];
            if (rhs_numbers.find(rhs) == rhs_numbers.end()) {
                rhs_numbers[rhs] = 0;
            }

            rhs_numbers[rhs] += 1;
        }

        for (auto pair : rhs_numbers) {
            out_file << "Number of b_i = " << pair.first << " : " << pair.second << std::endl;
        }
    }

    out_file.close();

    std::cout << "Generated statistics." << std::endl;
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
        for (auto perturbation : { "0.2" }) {
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
            solution.push_back(problem.columns_[j]->isInCurrentSolution() ? 1 : 0);
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
        characteristics << "Columns number : " << problem.columns_.size() << std::endl;
        characteristics << "Tasks number : " << problem.rhs_.size() << std::endl;
        characteristics << "Number of columns in solution : " << n_cols_in_solution << std::endl;
        characteristics << "Number of artificial variables in solution : " << n_ac_cols_in_solution << std::endl;
        characteristics << "Meannumber of non zeroes : " << ((double)n_non_zeros) / problem.columns_.size() << std::endl;
        for (auto pair : rhs_numbers) {
            characteristics << "Nombre de b_i = " << pair.first << " : " << pair.second << std::endl;
        }
        characteristics.close();
    }

}