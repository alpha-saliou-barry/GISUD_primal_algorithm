#pragma once

#include "ISUD_Base.h"
#include <ilcplex/ilocplex.h>

class IB_ReducedProblem {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;
	// Vector of activeConstraints in the RP
	std::vector <std::string> activeConstraints_;

public:
	// Constructor of reduced problem
	// "psolutionMethod" is a pointer on the problem
	// "activeConstraints" is the set of active constraints in the RP
	IB_ReducedProblem(ISUD_Base* psolutionMethod, std::vector <std::string> activeConstraints) {
		psolutionMethod_ = psolutionMethod;
		activeConstraints_ = activeConstraints;
	}
	

	// Solve RP and put the result in newSolution
	// "currentSolution" is the current ISUD solution
	// Returns new integral solution in "newSolution"
	// "previous_solution" contains previous solution cost
	// "max_size" is the max size of RP
	double solveProblem(std::vector<int>& currentSolution, std::vector<int>* newSolution,  double previous_solution, int phase = -1);


	// Get dual variables for the RP
	// Return dual variables in "duals"
	double getDuals(std::vector<double>* duals);
};