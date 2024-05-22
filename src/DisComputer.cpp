#include "DisComputer.h"
#include "IncompatibilityDegree.h"
#include "ComplementaryProblemDis.h"

#include <string>
std::map<int, std::vector<int>> DisComputer::computeDisaggregation(std::vector<int> colsIn, std::vector<int> colsOut, std::map<int, std::map<int, int>>* posColumnsAffectation, double* objective) {
	// On d�sagr�ge toutes les taches de la direction donn�e
	std::map<int, std::vector<int>> associations;
	std::set<int> involvedTasks;
	std::vector<int> involvedColumns;
	std::map<std::string, int> tasksCorrespondance;
	for (int i = 0; i < initial_problem->tasks_.size(); i++) {
		tasksCorrespondance[initial_problem->tasks_[i]] = i;
	}

	std::vector<IB_Column*> all_columns = initial_problem->columns_;
	int n_pos_columns = 0;
	for (int i = 0; i < all_columns.size(); i++) {
		if (all_columns[i]->isInCurrentSolution()) {
			involvedColumns.push_back(i);
			n_pos_columns += 1;
		}
	}

	std::cout << n_pos_columns << " colonnes initiales" << std::endl;

	for (int i = 0; i < colsIn.size(); i++) {
		if (!initial_problem->columns_[colsIn[i]]->isInCurrentSolution()) {
			std::cout << "colonne " << i << std::endl;
		}
		involvedColumns.push_back(colsIn[i]);
		for (auto contrib : all_columns[colsIn[i]]->getContribs()) {
			involvedTasks.insert(tasksCorrespondance[contrib.first]);
		}
	}

	for (int i = 0; i < colsOut.size(); i++) {
		for (auto contrib : all_columns[colsOut[i]]->getContribs()) {
			involvedTasks.insert(tasksCorrespondance[contrib.first]);
		}
	}

	for (int i = 0; i < initial_problem->tasks_.size(); i++) {
		associations[i] = std::vector<int>();
		if (involvedTasks.count(i)) {
			for (int j = 0; j < initial_problem->rhs_[i]; j++) {
				associations[i].push_back(1);
			}
		}
		else {
			associations[i].push_back(initial_problem->rhs_[i]);
		}
		
	}

	return associations;
}