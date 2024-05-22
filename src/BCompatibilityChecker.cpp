#include "BCompatibilityChecker.h"
#include <set>
#include <Eigen/Dense>

// Check if column is include in other column
bool BCompatibilityChecker::isColumnInclude(IB_Column* col1, IB_Column* col2) {
	std::set<std::string> cols2Contribs;
	for (auto pair: col2->getContribs()) {
		if (pair.second != 0) {
			cols2Contribs.insert(pair.first);
		}
	}

	for (auto pair : col1->getContribs()) {
		if (pair.second != 0 && !cols2Contribs.count(pair.first)) {
			return false;
		}
	}

	return true;
}

// Init the checking of binary compatible columns
void BCompatibilityChecker::init() {
	// Initialise les subset columns
	std::map<std::string, std::set<int>> tasksColumns;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		IB_Column* column = psolutionMethod_->columns_[i];
		for (auto contrib: column->getContribs()) {
			if (contrib.second != 0) {
				if (tasksColumns.find(contrib.first) == tasksColumns.end()) {
					tasksColumns[contrib.first] = std::set<int>();
				}

				tasksColumns[contrib.first].insert(i);
			}
		}

		subsetColumns_[i] = std::vector<int>();
	}

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		std::set<int> inColumns;
		for (auto contrib : psolutionMethod_->columns_[i]->getContribs()) {
			if (contrib.second != 0) {
				std::set<int> result;
				if (inColumns.size() == 0) {
					inColumns = tasksColumns[contrib.first];
				}

				std::set_intersection(inColumns.begin(), inColumns.end(), tasksColumns[contrib.first].begin(), tasksColumns[contrib.first].end(), std::inserter(result,
					result.begin()));
				inColumns = result;
			}
		}

		for (auto col : inColumns) {
			if (col != i) {
				subsetColumns_[col].push_back(i);
			}
		}
	}

	std::cout << "subset calcules" << std::endl;

	for (int i = 0; i < psolutionMethod_->columns_.size();i++) {
		IB_Column* column = psolutionMethod_->columns_[i];
		if (!column->isInCurrentSolution()) {
			updateCompatibilityStatus(i);
		}
	}
}

// Update compatibility status of column column_id

void BCompatibilityChecker::updateCompatibilityStatus(int column_id) {
	std::vector<int> indices;
	if (isBinaryCompatible(column_id, &indices)) {
		makeCompatibles_[column_id] = indices;
		psolutionMethod_->columns_[column_id]->setCompatible();
	}
	else {
		psolutionMethod_->columns_[column_id]->setIncompatible();
	}
}

// Check if the column column_id is binary compatible with columns in indices
bool BCompatibilityChecker::isBinaryCompatible(int column_id, std::vector<int>* indices) {
	indices->clear();
	
	std::vector<int> colsToCheck;
	for (auto pcolId : subsetColumns_[column_id]) {
		IB_Column* column = psolutionMethod_->columns_[pcolId];
		if (column->isInCurrentSolution() && std::find(colsToCheck.begin(), colsToCheck.end(), pcolId) == colsToCheck.end()) {
			colsToCheck.push_back(pcolId);
		}
	}
	 
	IB_Column* column = psolutionMethod_->columns_[column_id];
	std::set<std::string> tasksToCover;
	for (auto pair : column->getContribs()) {
		if (pair.second != 0) {
			tasksToCover.insert(pair.first);
		}
	}

	return checkCompatibility(colsToCheck, tasksToCover, indices);
}

// Check compatibility of columns colsToCheck. Do they cover the tasks tasksToCover ?
bool BCompatibilityChecker::checkCompatibility(std::vector<int> colsToCheck, std::set<std::string> tasksToCover, std::vector<int>* indices) {
	if (tasksToCover.size() == 0) {
		return true;
	}

	if (colsToCheck.size() == 0) {
		return false;
	}

	for (int i = 0; i < colsToCheck.size();i++) {
		int colId = colsToCheck[i];
		IB_Column* column = psolutionMethod_->columns_[colId];
		std::set<std::string> remainingTasks = (tasksToCover);

		bool inColumn = true;
		for (auto pair : column->getContribs()) {
			if (pair.second != 0) {
				if (!tasksToCover.count(pair.first)) {
					inColumn = false;
				}
				else {
					remainingTasks.erase(pair.first);
				}
			}
		}

		if (inColumn) {
			indices->push_back(colId);
			std::vector<int> remainingColsToCheck = colsToCheck;
			remainingColsToCheck.erase(remainingColsToCheck.begin(), remainingColsToCheck.begin() + i + 1);

			if (checkCompatibility(remainingColsToCheck, remainingTasks, indices)) {
				return true;
			}
			else {
				indices->pop_back();
			}
		}
	}

	return false;
}