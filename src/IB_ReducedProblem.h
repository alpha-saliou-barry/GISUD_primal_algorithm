#pragma once

#include "ISUD_Base.h"
#include <ilcplex/ilocplex.h>

class IB_ReducedProblem {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;
	// Vector of activeConstraints in the RP
	std::vector <std::string> activeConstraints_;
	// Optimality GAP of the RP (when to stop)
	double gap_;

public:
	// Constructor of reduced problem
	IB_ReducedProblem(ISUD_Base* psolutionMethod, std::vector <std::string> activeConstraints, double gap) {
		psolutionMethod_ = psolutionMethod;
		activeConstraints_ = activeConstraints;
		gap_ = gap;
	}
	
	// Solve RP and put the result in newSolution
	double solveProblem(std::vector<int>& currentSolution, std::vector<int>* newSolution,  double previous_solution, int phase = -1);

	// Get dual variables for the RP
	double getDuals(std::vector<double>* duals);
};