#pragma once
#include "IB_Column.h"
#include "ISUD_Base.h"
#include <set>

std::vector<std::pair<std::string, int>> read_tasks(std::string file);
std::vector<IB_Column*> read_columns(std::string file, double proba = -1);
std::set<std::string> read_initial_solution(std::string file);
ISUD_Base get_problem_random(std::string folder);
ISUD_Base get_problem(std::string folder);