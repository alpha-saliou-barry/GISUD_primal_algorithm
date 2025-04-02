#pragma once
#include <vector>
#include <string>
#include <fstream>
#include "Files.h"

// Get problem statistics of initial folders "initial_folders", with problem names : "problem_names"
void getStatistics(std::vector<std::string> initial_folders, std::vector<std::string> problem_names);