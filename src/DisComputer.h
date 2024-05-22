#pragma once
#include "DisProblem.h"

class DisComputer {
private:
	ISUD_Base* initial_problem;
	std::map<int, int> originalProblemColumns_;
	int phase_;
public:
	DisComputer(ISUD_Base* problem, std::map<int, int> originalProblemColumns, int phase) {
		initial_problem = problem;
		originalProblemColumns_ = originalProblemColumns;
		phase_ = phase;
	}

	std::map<int, std::vector<int>> computeDisaggregation(std::vector<int> colsIn, std::vector<int> colsOut, std::map<int, std::map<int, int>>* posColumnsAffectation, double* objective);
};