#pragma once
#include "ISUD_Base.h"
#include <set>
#include <string>
#include <algorithm>
#include "Files.h"
#include "ISUD.h"

// Perturbate the optimal solution of the problem problem at a given percentage
void perturbateFinalSolution(ISUD_Base* problem, double percentage = 0.5);

// Return artificial columns of cost "max_cost" to add to a problem on tasks "tasks" with a remaining rhs "remaining_rhs"
std::vector<IB_Column*> createInitialSolution(std::vector<std::string> tasks, std::vector<int> remaining_rhs, double max_cost);

// Generate problems from initial_folders to out_directories by perturbing the final solution
void generateProblems(std::vector<std::string> initial_folders, std::vector<std::string> out_directories);