#pragma once
#include <ilcplex/ilocplex.h>
#include "ISUD_Base.h"

class CplexMIP {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;

public:
	// Constructor
	CplexMIP(ISUD_Base* psolutionMethod) : psolutionMethod_(psolutionMethod) {

	}

	// Solve CPLEX MIP
	double solve(std::vector<int>& currentSolution, std::vector<int>* solution, bool relaxation = false, std::string path = "", std::string export_path = "");
	
	// Solve CPLEX MIP from path
	double solveFromFile(std::string path, std::string solution_file);
};