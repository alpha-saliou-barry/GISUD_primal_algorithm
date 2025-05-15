#pragma once

#include "ISUD_Base.h"
#include <vector>
#include <string>
#include <set>
#include <algorithm>
// Pertubate final solution of problem "problem" at a given percetage "percentage".
// "percentage" is the percentage of columns that remains in final solution.
void perturbateFinalSolution(ISUD_Base* problem, double percentage = 0.5);

// Add artificial columns
// "tasks" is the problem tasks
// "remaining_rhs" is the remaining rhs to cover
// "max_cost" is the new cost given to each columns
std::vector<IB_Column*> addArtificialColumns(std::vector<std::string> tasks, std::vector<int> remaining_rhs, double max_cost);
