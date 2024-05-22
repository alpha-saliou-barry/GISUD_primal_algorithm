#pragma once
#include "ISUD_Base.h"


class DisProblem {
private:
	ISUD_Base* psolutionMethod_;
	std::vector<int> initial_rhs;
	std::map<int, std::vector<int>> rhs_mapping;
	std::vector<int> involvedColumns_;
	std::vector<IB_Column*> allColumns;
	std::vector<int> initialColumns;
	std::map<std::pair<int, int>, std::string> tasksMappingInv;
	std::map<int, std::map<int, int>> columnsCorrespondance;
	std::map<int, std::vector<int>> initialColumnAssociatedInv;
	std::map<std::string, int> tasksCorrespondance;
	ISUD_Base* current_problem;
	

public:
	std::map<int, int> originalProblemColumns_;
	std::map<std::string, std::pair<int, int>> tasksMapping;
	std::map<int, int> initialColumnAssociated;
	std::vector<std::string> newTasks;
	std::vector<int> newRhs;
	DisProblem(ISUD_Base* problem,  std::vector<int> involvedColumns, std::map<int, std::vector<int>> mapping) {
		involvedColumns_ = involvedColumns;
		rhs_mapping = mapping;
		psolutionMethod_ = problem;
		for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
			tasksCorrespondance[psolutionMethod_->tasks_[i]] = i;
		}

	};

	ISUD_Base* getISUDBase();

	std::vector<std::map<int, int>> getTasksInstances(int indice, IB_Column* column, std::map<std::pair<int, int>, int>* remainingRhs, std::map<int, std::map<int, int>>* posColumnsAffectation = NULL);
	std::vector<std::map<int, int>> getTasksInstancesBis(std::vector<std::pair<std::string, int>> contribs);
	std::pair<int, int> getPosition(std::string taskName);
	std::vector<int> getAssociatedColumns(int task, int n, int column_no);
	std::vector<IB_Column*> getAllColumns() {
		return allColumns;
	}
	int constructColumns(std::map<int, int> originalProblemColumns, int phase = 1000, std::map<int, std::map<int, int>>* posColumnsAffectation = NULL);
	std::vector<double> aggregate(std::vector<double>& solution, int n_cols);
	void deleteProblem();
};