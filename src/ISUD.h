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

	// Return if the support of one direction "solution" can be included in an integer direction.
	// "column" correspond to the artificial column and "acolId" to the id of the artificial column.
	bool canBeInIntegerDirection(IB_Column* newColumn, std::vector<double>& solution, std::vector<int>* support = NULL, int acolId = -1);

	// Return if the solution "solution" is integral 
	static bool isIntegral(std::vector<double>& cpSolution);

	// Constructor of GISUD
	// parameters are the gap "gap_", a boolean to enable the column addition strategy "addColumns", a boolean to know if we should check the binary compatibility of columns "checkBinaryCompatibility" and a boolean to know if we want to compare column addition strategy and zoom strategy "compete_"
	ISUD(ISUD_Base* problem, double gap_, bool addColumns = true, bool checkBinaryCompatibility = true, bool compete_ = false);
	
	// Return a binary compatible column of negative reduced cost in "colsIn", "colsOut"
	bool getBCompatibleColumn(std::vector<int>* colsIn, std::vector<int>* colsOut);
	
	// Return a binary compatible column of negative reduced cost with incompatibility degree in "colsIn", "colsOut"
	bool getBCompatibleColumnId(std::vector<int>* colsIn, std::vector<int>* colsOut);
	
	// Solve GSPP in folder : "path"
	void solve(std::string path = "");
	
	// Pivot columns "colsIn", "colsOut" in the solution, recompute compatibilities if "recomputeCompatibilities" is true
	void pivotColumnsInSolution(std::vector<int>& colsIn, std::vector<int>& colsOut, bool recompute = true);
	
	// ZOOM procedure from phase "isudPhase", with initial zoom solution "solution". It returns the solution in "colsIn" and "colsOut".
	bool zoom(int, std::vector<double>& solution, std::vector<int>* colsIn, std::vector<int>* colsOut);
	
	// Complementary problem in phase "initial_phase" with column addition.
	// The artificial column id is in "acolId", aggregating columns in "colsIn" "colsOut"
	// Previous objective in "previous_objective" and number of artificial columns created in "n_a_cols"
	// Put solution in nColsIn, nColsOut
	// Returns a boolean (if the function found an integral direction or not) and an integer (number of artificial columns created).

	std::pair<bool, int> cpWithArtificialColumn(int acolId, std::set<int> colsIn, std::set<int> colsOut,
		std::vector<int>& ncolsIn, std::vector<int>& ncolsOut, int initial_phase, int n_a_cols = 1, double previous_objective = 0);
	
	// Return sequence of GISUD phase (incompatibility degrees)
	std::vector<int> getPhaseSequence();

	// Add a line to the output
	void addLine(double newCost, double amelioration, int n_iterations, int n_added_columns, int n_success_add_columns, bool has_zoom, int phaseMax, bool addedColumn, int n_added_columns_it,
		int pivot_distance, double iteration_time, double global_time);
	// Search for an integer direction in the support of one direction "in_columns".
	// Return the integer direction in "colsIn", "colsOut"
	bool searchSubDirection(std::vector<int> in_columns, std::vector<int>* out_columns, std::vector<int>* colsIn);
	
	// Add row for the comparison of ZOOM and column addition
	void addCompeteRow(double amelioration_rc, int time_rc, double amelioration_zoom, int time_zoom, int disp_size = 0, double last_objective = 0, double remaining = 0); 
};

// Compute incompatibility degrees of the columns 'columns'
void calcIncompatibilityDegrees(IncompatibilityDegree id, std::vector<IB_Column*> columns);

// Compute incompatibility degrees of the columns "columns_to_recompute" with "n_threads" threads
void computeIncompatibilityDegreeWithThreads(std::vector<IB_Column*> columns_to_recompute, ISUD_Base* psolutionMethod_, int n_threads = 8);