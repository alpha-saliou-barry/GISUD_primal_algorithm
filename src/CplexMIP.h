#pragma once
#include <ilcplex/ilocplex.h>
#include "ISUD_Base.h"

class CplexMIP {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;

public:
	// Constructor of CplexMIP
	// Pointer of the problem in "psolutionMethod"
	CplexMIP(ISUD_Base* psolutionMethod) : psolutionMethod_(psolutionMethod) {

	}

	// Solve GSPP with CPLEX
    // Gap is the optimality gap
	// Initial solution in "currentSolution"
	// Put solution in "solution"
	// "Relaxation" is a boolean to solve only the linear relaxation
	// "path" is the path of CPLEX output (optional)
	// "export_path" is to export CPLEX model (optional)

	double solve(std::vector<int>& currentSolution, std::vector<int>* solution, bool relaxation = false, double gap = -1, 
		std::string path = "", std::string export_path = "");
};