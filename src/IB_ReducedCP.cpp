#include "IB_ReducedCP.h"
#include <ilcplex/ilocplex.h>
#include <cmath>
#include "IB_CompatibilityChecker.h"


// Solve reduced CP for ZOOM
// Current objective of GISUD in "currentValue"
// Dual variables of the reduced CP in "dualVariables"
// Active constraints of the reduced problem in "activeConstraintsRP"
// Put columns to add to the Reduced Problem in "inColumns"
// Put columns indices of "inColumns" in "inColumnsIndices"
// "n_calls" is the number of calls to the CP problem
// Maximum incompatibility degree of columns of the reduced CP in "phase"
bool IB_ReducedCP::solve(double currentValue, std::vector<double> dualVariables, std::vector<std::string> activeConstraintsRP, std::vector<IB_Column*>* inColumns, std::vector<int>* inColumnsIndices, int n_calls, int phase) {
	// Problem construction
	IloEnv env;
	IloModel mod(env);
	IloNumVarArray vars(env);
	IloObjective cost = IloAdd(mod, IloMinimize(env));
	std::map < std::string, int > constraintsIds;
	std::vector<int> colsIndices;

	

	IloRangeArray constraints(env, psolutionMethod_->tasks_.size() + 1, 0, 0);
	constraints[psolutionMethod_->tasks_.size()].setBounds(1, 1);
	
	for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
		constraintsIds[psolutionMethod_->tasks_[i]] = i;
	}
	
	for (int i = 0; i < constraints.getSize(); i++) {
		mod.add(constraints[i]);
	}
	
	std::vector<double> completeDualVariables;
	for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
		completeDualVariables.push_back(0);
	}
	for (int i = 0; i < dualVariables.size(); i++) {
		completeDualVariables[constraintsIds[activeConstraintsRP[i]]] = dualVariables[i];
	}
	
	bool oneIncompatibleColumn = false;

	

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {

		IB_Column* column = psolutionMethod_->columns_[i];

		if ((!column->isInP() && column->getState() == IB_Column::STATE_INCOMPATIBLE && (phase == -1 || column->getPhase() <= phase)) || column->isInPIndependent()) {
			int sign = column->isInPIndependent() ? -1 : 1; 

			double cost_ = sign * column->getCost();
			if (column->isInPIndependent()) {
				cost_ = 0;
				for (auto contrib : column->getContribs()) {
					if (contrib.second != 0) {
						cost_ +=  sign * contrib.second * completeDualVariables[constraintsIds[contrib.first]];
					}

				}
			}
			
			IloNumColumn col = cost(cost_);
			for (auto contrib : column->getContribs()) {
				if (contrib.second != 0) {
					col += constraints[constraintsIds[contrib.first]](sign* contrib.second);
				}

			}

			

			if (!column->isInP() && column->getState() == IB_Column::STATE_INCOMPATIBLE) {
				oneIncompatibleColumn = true;
			}


			if (column->isInPIndependent()) {
				vars.add(IloNumVar(col, -1000000, IloInfinity));
			}
			else {
				col += constraints[psolutionMethod_->tasks_.size()](1);
				vars.add(IloNumVar(col, 0, IloInfinity));
			}
			

			colsIndices.push_back(i);

			col.end();
		}
	}



	std::cout << vars.getSize() << " variables in phase " << phase << std::endl;
	IloCplex cplex(mod);
	cplex.setParam(IloCplex::Param::Simplex::Display, 0);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-2);
	cplex.setParam(IloCplex::Param::Simplex::Limits::LowerObj, -1);

	for (int i = 0; i < n_calls; i++) {
		IloNumArray vals(env);
		bool success = cplex.solve();
		if (!success || cplex.getObjValue() >= 0 || fabs(cplex.getObjValue() * psolutionMethod_->sum_bi) / currentValue <= 0.1) { 
			// std::cout << "The success is : " <<  (success ? cplex.getObjValue() : -1) << std::endl;
			return (inColumns->size() != 0);
		}

		cplex.getValues(vars, vals);
		double sum = 0;
		std::set<std::string> coveredTasks;
		for (int j = 0; j < vals.getSize(); j++) {
			if (j < vals.getSize()) {
				if (fabs(vals[j]) >= 1e-5) {

					if (!psolutionMethod_->columns_[colsIndices[j]]->isInP() && psolutionMethod_->columns_[colsIndices[j]]->getState() == IB_Column::STATE_INCOMPATIBLE) {
						inColumns->push_back(psolutionMethod_->columns_[colsIndices[j]]);
						inColumnsIndices->push_back(colsIndices[j]);
					}
					else {
						sum += vals[j];
					}

					for (auto contrib : psolutionMethod_->columns_[colsIndices[i]]->getContribs()) {
						coveredTasks.insert(contrib.first);
					}
				}

			}
			else {
				std::cout << "Value of artificial variable : " << vals[j] << std::endl;
			}
		}

		int n2 = 0;
		// Deletion of columns
		for (int k = 0; k < vars.getSize(); k++) {
			bool delete_ = false;
			for (auto task : coveredTasks) {
				if (psolutionMethod_->columns_[colsIndices[k]]->findContribution(task) > 1e-4) {
					delete_ = true;
				}
			}

			if (fabs(vals[k]) > 1e-5) {
				vars[k].setBounds(0, 0);
				n2 += 1;
			}
		}

		std::cout << n2 << " variables put to 0" << std::endl;
		// Deletion of constraints
		for (auto task : coveredTasks) {
			//mod.remove(constraints[constraintsIds[task]]);
		}

		std::cout << sum << std::endl;
	}
		env.end();


		return true;
	}
