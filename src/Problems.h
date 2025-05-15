#pragma once
#include "ISUD_Base.h"
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include "InitialSolutions.h"
#include "CplexMIP.h"
#include "ISUD.h"

// Generate all instances
void generateProblems();


// Construct ISUD_Base from paths
// "mps_file" is the mps_file path
// Compute initial solution for problem "problem"
// "solution" is the optimal solution of problem "problem"
// "perturbation_percent" is the percentage of remaining columns in new initial solution
// Out dir is the output directory where to put new problems with new initial solutions
int computeInitialSolution(std::vector<int> solution, double perturbation_percent, ISUD_Base& problem, std::string out_dir);

// Export problem to a file in the GISUD format
// "problem" is the problem
// "columns_file" is the path to the columns file
// "rhs_file" is the path to the rhs file
// "initial_solution_file" is the path of the initial solution file
// "fixed_cost_path" is the path of the fixed cost file
void export_problem(ISUD_Base* problem, std::string columns_file, std::string rhs_file, std::string initial_solution_file, std::string fixed_cost_file = "");

// Construct ISUD_Base from paths
// "mps_file" is the mps_file path
// "rhs_file" is the rhs_file path
// "initial_solution_file" is the initial solution path
// "fixed_cost_file" is the fixed cost file
ISUD_Base generateProblemFromMps(std::string mps_file, std::string rhs_file, std::string initial_solution_file = "", std::string fixed_cost_file = "");