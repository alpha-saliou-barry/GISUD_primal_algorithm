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
	// Constructor of reduced problem. "psolutionMethod" is a pointer on the problem and "activeConstraints" is the vector of active constraints in RP
	IB_ReducedProblem(ISUD_Base* psolutionMethod, std::vector <std::string> activeConstraints) {
		psolutionMethod_ = psolutionMethod;
		activeConstraints_ = activeConstraints;
	}

	// Solve RP from "currentSolution" and put the result in "newSolution". The previous objective is "previous_solution". "phase" is the phase of the RP.
	double solveProblem(std::vector<int>& currentSolution, std::vector<int>* newSolution,  double previous_solution, int phase = -1);
	
	// Get dual variables for the RP. Put the result in "duals".
	double getDuals(std::vector<double>* duals);
};