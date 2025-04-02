#include "Statistics.h"


// Get problem statistics of initial folders "initial_folders", with problem names : "problem_names"
void getStatistics(std::vector<std::string> initial_folders, std::vector<std::string> problem_names, std::string out_path) {
    std::ofstream out_file;
    out_file.open(out_path);

    // Get folder statistics
    for (int k = 0; k < initial_folders.size(); k++) {
        std::string initial_folder = initial_folders[k];
        ISUD_Base problem = constructISUDProblem(initial_folder + "/columns.txt",
            initial_folder + "/rhs.txt", initial_folder + "/initial.txt");

        int n_non_zeros = 0;
        for (auto column : problem.columns_) {
            n_non_zeros += column->getContribs().size();
        }

        out_file << std::endl << problem_names[k] << std::endl;
        out_file << "Initial rows number: " << problem.tasks_.size() << std::endl;
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

    std::cout << "Statistics generated" << std::endl;
}