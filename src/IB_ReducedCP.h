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
	IB_ReducedCP(ISUD_Base* psolutionMethod, std::set<int> activePColumns, std::vector<std::string> activeConstraints, std::vector<int> columnsIndices, std::vector<Eigen::VectorXf*> eigenVectors,
		std::vector<double> costs) : psolutionMethod_(psolutionMethod), columnsIndices_(columnsIndices), eigenVectors_(eigenVectors),
		costs_(costs), activeConstraints_(activeConstraints), activePColumns_(activePColumns) {

	}

	// Solve the reducedCP
	bool solve(double currentValue, std::vector<double> dualVariables, std::vector<std::string> activeConstraintsRP, std::vector<IB_Column*>* inColumns, std::vector<int>* inColumnsIndices, int n_calls=5, int phase=-1);
};