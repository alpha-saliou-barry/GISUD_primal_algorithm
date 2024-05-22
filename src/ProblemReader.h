#pragma once
#include "IB_Column.h"
#include "ISUD_Base.h"
#include <set>

// Read tasks from file
std::vector<std::pair<std::string, int>> read_tasks(std::string file);
// Read columns from file
std::vector<IB_Column*> read_columns(std::string file, double proba = -1);
// Read initial solution from file
std::set<std::string> read_initial_solution(std::string file);
// Get random problem
ISUD_Base get_problem_random(std::string folder);
// Get problem from folder
ISUD_Base get_problem(std::string folder);