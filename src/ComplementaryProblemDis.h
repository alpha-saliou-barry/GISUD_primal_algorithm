#pragma once

#include "ISUD_Base.h"

#include <ilcplex/ilocplex.h>


class ComplementaryProblemDis {
private:
	ISUD_Base* psolutionMethod_;
	std::vector<double> normalizationConstraint_;
	int current_phase;
	IloEnv env;
	IloModel mod;
	IloCplex cplex;
	IloNumVarArray vars;
	std::map<int, int> varsIds;
	IloRangeArray constraints;
	std::map<std::string, int> constraintsIds;
	std::vector<int> colsIndices;
	IloObjective cost;
public:
	ComplementaryProblemDis(ISUD_Base* psolutionMethod, int phase = -1) : psolutionMethod_(psolutionMethod),
		mod(env), vars(env), constraints(env, psolutionMethod_->tasks_.size() + 1, 0, 0), current_phase(phase), cplex(mod)
	{
		for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
                  normalizationConstraint_.push_back(psolutionMethod_->columns_[i]->isInCurrentSolution() ? 1:(psolutionMethod_->columns_[i]->getPhase() + 1));
                  //normalizationConstraint_.push_back(1);
		}
		cost = IloAdd(mod, IloMinimize(env));
	}

	void constructProblem();

	double solve(std::vector<double>* solution, std::map<int, int> originalProblemColumns, 
	std::map<std::string, std::pair<int, int>> tasksMapping, std::vector<double>& duals, bool first = false, bool export_ = false);
	void destroy();
};
