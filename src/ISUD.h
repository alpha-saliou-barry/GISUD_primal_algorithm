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
	// If we want to use disaggregation
	bool disEnabled;
	// If we want to compare strategy (disaggregation, column addition) with zoom
	bool compete;
	// Size of the disaggregated problem
	std::string dis_problem_size = "";
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

	//Return if the support of one direction can be included in an integer direction
	bool canBeInIntegerDirection(IB_Column* newColumn, std::vector<double>& solution, std::vector<int>* support = NULL, int acolId = -1);

	// Return if the solution is integral
	static bool isIntegral(std::vector<double>& cpSolution);

	// Constructor of ISUD, addColumns strategy or disaggregationStrategy. Compete true means we want to compare the strategy (addColumns or disaggregation) with zoom.
	ISUD(ISUD_Base* problem, bool addColumns = true, bool checkBinaryCompatibility = true, bool disE= false, bool compete_ = false);

	// Return binary compatible column of negative reduced cost
	bool getBCompatibleColumn(std::vector<int>* colsIn, std::vector<int>* colsOut);

	// Return binary compatible column of negative reduced cost with incompatibility degree
	bool getBCompatibleColumnId(std::vector<int>* colsIn, std::vector<int>* colsOut);

	// Main procedure of ISUD, solve the problem and stock the output to path
	void solve(std::string path = "");

	// Pivot columns colsIn, colsOut in the solution, recompute compatibilities if recomputeCompatibilities is true
	void pivotColumnsInSolution(std::vector<int>& colsIn, std::vector<int>& colsOut, bool recompute = true);

	// ZOOM procedure
	bool zoom(int, std::vector<int> seqPhases, std::vector<double>& solution, std::vector<int>* colsIn, std::vector<int>* colsOut);

	// Complementary problem with column addition
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