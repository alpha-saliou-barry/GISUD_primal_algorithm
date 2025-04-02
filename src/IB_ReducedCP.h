#pragma once

#include "ISUD_Base.h"
#include <Eigen/Dense>

class IB_ReducedCP {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;
	// Active constraints in the CP
	std::vector<std::string> activeConstraints_;
	// Vector of original columns indices of the CP
	std::vector<int> columnsIndices_;
	// Costs of columns in the CP
	std::vector<double> costs_;
	// Positive columns indices in the CP
	std::set<int> activePColumns_;

public:
	// Constructor of ReducedCP.
	// Pointer of the problem in "psolutionMethod_"
	// Active positive columns in "activePColumns"
	// Active CP constraints in "activeConstraints"
	// Columns indices of the CP in "columnsIndices"
	IB_ReducedCP(ISUD_Base* psolutionMethod, std::set<int> activePColumns, std::vector<std::string> activeConstraints, std::vector<int> columnsIndices, std::vector<Eigen::VectorXf*> eigenVectors,
		std::vector<double> costs) : psolutionMethod_(psolutionMethod), columnsIndices_(columnsIndices),
		costs_(costs), activeConstraints_(activeConstraints), activePColumns_(activePColumns) {

	}

	// Solve reduced CP for ZOOM
	// Current objective of GISUD in "currentValue"
	// Dual variables of the reduced CP in "dualVariables"
	// Active constraints of the reduced problem in "activeConstraintsRP"
	// Put columns to add to the Reduced Problem in "inColumns"
	// Put columns indices of "inColumns" in "inColumnsIndices"
	// "n_calls" is the number of calls to the CP problem
	// Maximum incompatibility degree of columns of the reduced CP in "phase"

	bool solve(double currentValue, std::vector<double> dualVariables, std::vector<std::string> activeConstraintsRP, std::vector<IB_Column*>* inColumns, std::vector<int>* inColumnsIndices, int n_calls=5, int phase=-1);
};