#pragma once
#include "ISUD_Base.h"
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include "InitialSolutions.h"

// Generate ISUD_Base problem from GISUD format
// Path of columns in "columns_path"
// Path of right-hand sides in "rhs_path"
// Path of initial solution in "initial_solution_path"
// Path of fixed cost file in "fixed_cost_path"
ISUD_Base constructISUDProblem(std::string columns_path, std::string rhs_path, std::string initial_solution_path, std::string fixed_cost_path = "");


// Export ISUD_Base problem to GISUD format
// Put columns in "columns_path" file
// Put right-hand-sides in "rhs_path" file
// Put fixed cost in "fixed_cost_path"
void export_problem(ISUD_Base* problem, std::string columns_path, std::string rhs_path, std::string initial_solution_path, std::string fixed_cost_path = "");
