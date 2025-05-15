#pragma once

#include "ISUD_Base.h"
#include "BCompatibilityChecker.h"
#include <iostream>
#include <fstream>
#include "IncompatibilityDegree.h"

class ISUD {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;
	// Current solution cost
	double currentCost_;
	// Current solution gap
	double current_gap;
	// Best lower bound
	double bestBound;
	// If we want to compare strategy (disaggregation, column addition) with zoom
	bool compete;
	// Binary compatibility checker
	BCompatibilityChecker bcompatibilityChecker_;
	// Indices of columns in currentSolution
	std::vector<int> currentSolution_;
	// Identifier of artificial column
	int artificial_column_id = 0;
	// Add columns strategy and if we want to check the binary compatibility exactly (without incompatibility degree)
	bool addColumns_, checkBinaryCompatibility_;
	// The file of the output (without or with competition)
	std::ofstream file, fileC;
	// Separator of the output
	const char separator = ' ';
	// Width of the column name in the output
	const int nameWidth = 30;
	// With of column value in the output
	const int numWidth = 10;
	// Gap when to stop
	double gap = 0.01;

public:
	// Return current solution
	std::vector<int> getCurrentSolution() {
		return currentSolution_;
	}

	// Return true if the support of one direction can be included in an integer direction
	// "solution" is the direction
	// "support" is the support of the artificial column
	// "acolId" is the artificial column id
	// New artificial columnis returned in column "column"
	bool canBeInIntegerDirection(IB_Column* column, std::vector<double>& solution, std::vector<int>* support = NULL, int acolId = -1);

	// Returns true if the solution "solution" is integral
	static bool isIntegral(std::vector<double>& solution);

	// Constructor of ISUD
	// "problem" is the pointer on the problem
	// "addColumns" is a boolean that is true if we want to enable column addition strategy
	// "checkBinaryCompatibility" is a boolean that is true if we want to check the exact binary compatibility of columns
	// "compete" is a boolean that is true if we want to compare column addition strategy and ZOOM

	ISUD(ISUD_Base* problem, bool addColumns, bool checkBinaryCompatibility, bool compete_);

	// Return binary compatible column of negative reduced cost
	// Return "colsIn" and "colsOut"
	// "colsIn" contains the binary compatible column id
	// "colsOut" contains the columns to remove from P
	bool getBCompatibleColumn(std::vector<int>* colsIn, std::vector<int>* colsOut);

	// Return binary compatible column of negative reduced cost by incompatibility degree
	// Return "colsIn" and "colsOut"
	// "colsIn" contains the binary compatible column id
	// "colsOut" contains the columns to remove from P
	bool getBCompatibleColumnId(std::vector<int>* colsIn, std::vector<int>* colsOut);

	// Main procedure of ISUD, solve the problem and stock the output to path
	void solve(std::string path = "");

	// Pivot columns "colsIn", "colsOut" in the solution, recompute compatibilities if recomputeCompatibilities is true
	void pivotColumnsInSolution(std::vector<int>& colsIn, std::vector<int>& colsOut, bool recompute = true);

	// ZOOM procedure
	// ZOOM phase is "isudPhase"
	// "solution" if the CP solution
	// Returns integral direction in "colsIn" (columns to remove from P) and "colsOut" (columns to enter in P)
	bool zoom(int, std::vector<double>& solution, std::vector<int>* colsIn, std::vector<int>* colsOut);

	// Complementary problem with column addition strategy
	// "acolId" is the artificial column id
	// "colsIn" and "colsOut" is the support of the artificial columns
	// Returns integral direction in "ncolsIn", "ncolsOut"
	// "initial_phase" is the phase of the CP
	// "n_a_cols" is the number of artificial columns added
	std::pair<bool, int> cpWithArtificialColumn(int acolId, std::set<int> colsIn, std::set<int> colsOut,
		std::vector<int>& ncolsIn, std::vector<int>& ncolsOut, int initial_phase, int n_a_cols = 1, double previous_objective = 0);

	// Return phase sequence for multiphase strategy
	std::vector<int> getPhaseSequence();

	// Add a line to the output
	void addLine(double newCost, double amelioration, int n_iterations, int n_added_columns, int n_success_add_columns, bool has_zoom, int phaseMax, bool addedColumn, int n_added_columns_it,
		int pivot_distance, double iteration_time, double global_time);

	// Search for an integer direction in the support of one direction
	bool searchSubDirection(std::vector<int> in_columns, std::vector<int>* out_columns, std::vector<int>* colsIn);

	// Add row for the comparison of ZOOM and column addition
	void addCompeteRow(double amelioration_rc, int time_rc, double amelioration_zoom, int time_zoom, int disp_size = 0, double last_objective = 0, double remaining = 0);  

	// Complementary problem with task disaggregation
	std::pair<bool, int> cpWithDisaggregation(std::vector<double>& duals, int phase, std::vector<int>* colsIn, std::vector<int>* colsOut, std::vector<int>& acolsOut, std::vector<int>& acolsIn, ISUD_Base* problem, std::map<int, int> originalProblemColumns, double past_objective, int n_cols,
	int* size_dis_problem, double bound);
};

// Compute incompatibility degrees of the columns 'columns'
void calcIncompatibilityDegrees(IncompatibilityDegree id, std::vector<IB_Column*> columns);