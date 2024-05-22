#include "IB_ReducedProblem.h"
#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

const int MAX_RP_SIZE = 25000;

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T>& v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	std::stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

// Solve RP and put the result in newSolution
double IB_ReducedProblem::solveProblem(std::vector<int>& currentSolution, std::vector<int>* newSolution, double previous_solution, int max_size) {
	//Construction du probl�me
	IloEnv env;
	IloModel mod(env);
	IloNumVarArray vars(env);

	//Valeurs initiales
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);

	std::vector <int> constraintsIds;
	int current_constraint = 0;

	for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
		if (current_constraint < activeConstraints_.size() && psolutionMethod_->tasks_[i] == activeConstraints_[current_constraint]) {
			current_constraint += 1;
			constraintsIds.push_back(i);
		}
	}

	std::map<std::string, int> tasksIndexes;
	for (int i = 0; i < activeConstraints_.size(); i++) {
		tasksIndexes[activeConstraints_[i]] = i;
	}



	// D�finition des contraintes
	IloRangeArray constraints(env, activeConstraints_.size(),  0, 0);
	std::cout << activeConstraints_.size() << " " << constraintsIds.size() << " contraintes" << std::endl;

	for (int i = 0; i < activeConstraints_.size(); i++) {
		constraints[i].setBounds(psolutionMethod_->rhs_[constraintsIds[i]], psolutionMethod_->rhs_[constraintsIds[i]]);
		mod.add(constraints);
	}

	IloObjective cost = IloAdd(mod, IloMinimize(env));
	cost.setConstant(psolutionMethod_->fixed_cost_);
	
	std::vector<int> varsIds;

	// Rajout des variables
	std::vector<int> phases;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		phases.push_back(!psolutionMethod_->columns_[i]->isInP() ? psolutionMethod_->columns_[i]->getPhase() : -1);
	}
	int n_columns = 0;

	for (auto i: sort_indexes(phases)) {
		IB_Column* column = psolutionMethod_->columns_[i];
		if (column->isInP() || column->getState() == IB_Column::STATE_COMPATIBLE_LINEARLY) {
			IloNumColumn col = cost(column->getCost());
			for (auto contrib : column->getContribs()) {
				if (tasksIndexes.find(contrib.first) != tasksIndexes.end()) {
					col += constraints[tasksIndexes[contrib.first]](contrib.second);
				}
			}

			vars.add(IloBoolVar(col));
			startVar.add(vars[vars.getSize() - 1]);
			startVal.add((currentSolution)[i]);
			varsIds.push_back(i);

			col.end();
			n_columns += 1;

			if (n_columns >= max_size) {
				break;
			}
		}
	}
	
	std::cout << vars.getSize() << " variables" << std::endl;
	// R�solution du probl�me
	IloCplex cplex(mod);
	cplex.setParam(IloCplex::Param::Simplex::Display, 0);
	cplex.setParam(IloCplex::Param::Threads, 8);
	cplex.setOut(env.getNullStream());
	cplex.addMIPStart(startVar, startVal);
	//cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.01);
	cplex.setParam(IloCplex::Param::TimeLimit, 3600 * 5);
	cplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 1);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);

	double objective;

	bool success = cplex.solve();
	int current_number_of_solutions = 0;
	double po = 2*previous_solution;
	while (success) {
		objective = cplex.getObjValue();
		int nb_sols = cplex.getSolnPoolNsolns();
		std::cout << nb_sols << " solutions." << std::endl;
		if (nb_sols == current_number_of_solutions || nb_sols >= 15) {
			IloNumArray vals(env);
			cplex.getValues(vars, vals, 0);
			for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
				newSolution->push_back(0);
			}

			for (int i = 0; i < vals.getSize(); i++) {
				newSolution->at(varsIds[i]) = vals[i] > 1e-4 ? 1 : 0;;
			}

			success = false;
		}
		else {
			current_number_of_solutions = nb_sols;
			cplex.solve();
			po = objective;
		}
	}

	env.end();

	return objective;
}

// Get dual variables for the RP
double IB_ReducedProblem::getDuals(std::vector<double>* duals) {
	//Construction du probl�me
	IloEnv env;
	IloModel mod(env);
	IloNumVarArray vars(env);

	//Valeurs initiales
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);

	std::vector <int> constraintsIds;
	int current_constraint = 0;

	for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
		if (current_constraint < activeConstraints_.size() && psolutionMethod_->tasks_[i] == activeConstraints_[current_constraint]) {
			current_constraint += 1;
			constraintsIds.push_back(i);
		}
	}



	// D�finition des contraintes
	IloRangeArray constraints(env, activeConstraints_.size(), 0, 0);
	std::cout << activeConstraints_.size() << " " << constraintsIds.size() << " contraintes" << std::endl;

	for (int i = 0; i < activeConstraints_.size(); i++) {
		constraints[i].setBounds(psolutionMethod_->rhs_[constraintsIds[i]], psolutionMethod_->rhs_[constraintsIds[i]]);
		mod.add(constraints);
	}

	IloObjective cost = IloAdd(mod, IloMinimize(env));

	std::vector<int> varsIds;

	// Rajout des variables
	std::vector<int> phases;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		phases.push_back(!psolutionMethod_->columns_[i]->isInP() ? psolutionMethod_->columns_[i]->getPhase() : -1);
	}
	int n_columns = 0;

	for (auto i : sort_indexes(phases)) {
		IB_Column* column = psolutionMethod_->columns_[i];
		if (column->isInP() || column->getState() == IB_Column::STATE_COMPATIBLE_LINEARLY) {
			IloNumColumn col = cost(column->getCost());
			for (int j = 0; j < activeConstraints_.size(); j++) {
				int contrib = column->findContribution(activeConstraints_[j]);
				if (contrib != 0) {
					col += constraints[j](contrib);
				}
			}

			vars.add(IloNumVar(col, 0, 1));
			varsIds.push_back(i);

			col.end();
			n_columns += 1;
		}
	}

	std::cout << vars.getSize() << " variables" << std::endl;
	// R�solution du probl�me
	IloCplex cplex(mod);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Param::TimeLimit, 3600 * 5);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);

	double objective = 0;

	bool success = cplex.solve();
	int current_number_of_solutions = 0;
	if(success) {
		objective = cplex.getObjValue();
		IloNumArray vals(env);			
		cplex.getDuals(vals, constraints);
		for (int i = 0; i < vals.getSize(); i++) {
			duals->push_back(vals[i]);
		}
	}

	env.end();

	std::cout << objective << " est l'objectif" << std::endl;

	return objective;
}
