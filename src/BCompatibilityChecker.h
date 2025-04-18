#pragma once
#include "ISUD_Base.h"
#include <set>

class BCompatibilityChecker {
private:
	// Pointer on the problem
	ISUD_Base* psolutionMethod_;



public:

	// Subset columns
	std::map<int, std::vector<int>> subsetColumns_;
	// Subsets that make the column compatibles
	std::map<int, std::vector<int>> makeCompatibles_;

	// Constructor of BCompatibilityChecker.
	// "pSolutionMethod" is a pointer on the problem
	BCompatibilityChecker(ISUD_Base* psolutionMethod) : psolutionMethod_(psolutionMethod) {

	}

	// Init the checking of binary compatible columns
	void init();

	// Update compatibility status of column "column_id"
	void updateCompatibilityStatus(int column_id);

	// Check if the column "column_id" is binary compatible with columns in "indices"
	bool isBinaryCompatible(int column_id, std::vector<int>* indices);

	// Check compatibility of columns "colsToCheck". Do they cover the tasks "tasksToCover" ?
	bool checkCompatibility(std::vector<int> colsToCheck, std::set<std::string> tasksToCover, std::vector<int>* indices);

	// Is "col1" included in "col2" ? (by tasks)
	static bool isColumnInclude(IB_Column* col1, IB_Column* col2);
};