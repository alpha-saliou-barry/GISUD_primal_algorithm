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
	// Vectors of columns in the CP
	std::vector<Eigen::VectorXf*> eigenVectors_;
	// Costs of columns in the CP
	std::vector<double> costs_;
	// Columns positive indices in the CP
	std::set<int> activePColumns_;

public:
	//Constructor of reducedCP
	// "psolutionMethod" is a pointer on the problem
	// "activePColumns" is the set of active columns of P
	// "activeConstraints"is the set of active constraints of ReducedCP
	// "columnsIndices" is the set of reduced CP columns indices
	// "eigenVectors" is the reduced CP columns
	// "costs" is the cost of the reducedCP columns
	IB_ReducedCP(ISUD_Base* psolutionMethod, std::set<int> activePColumns, std::vector<std::string> activeConstraints, std::vector<int> columnsIndices, std::vector<Eigen::VectorXf*> eigenVectors,
		std::vector<double> costs) : psolutionMethod_(psolutionMethod), columnsIndices_(columnsIndices), eigenVectors_(eigenVectors),
		costs_(costs), activeConstraints_(activeConstraints), activePColumns_(activePColumns) {

	}


	// Solve the reducedCP
	// "currentValue" is the ISUD current solution cost
	// "dualVariables" is the dual variables of problem RP
	// "activeConstraintsRP" is the set of RP active constraints
	// Return columns to add to the RP in "inColumns". Their column indices are in "inColumnsIndices"
	// "n_calls" is the number of calls to CP
	// "phase" is the phase of the CP
	bool solve(double currentValue, std::vector<double> dualVariables, std::vector<std::string> activeConstraintsRP, std::vector<IB_Column*>* inColumns, std::vector<int>* inColumnsIndices, int n_calls=5, int phase=-1);
};