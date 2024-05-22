#include "ComplementaryProblemDis.h"

void ComplementaryProblemDis::constructProblem() {
	int n_tasks = psolutionMethod_->tasks_.size();
	constraints[n_tasks].setBounds(1, 1);
	for (int i = 0; i < n_tasks + 1; i++) {
		mod.add(constraints[i]);
		if (i < psolutionMethod_->tasks_.size()) {
			constraintsIds[psolutionMethod_->tasks_[i]] = i;
		}
	} 

	std::map<int, int> constraintsIndexesAC;
	
        
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		IB_Column* column = psolutionMethod_->columns_[i];

		if (current_phase == -1 || (column->isInCurrentSolution() || column->getPhase() <= current_phase)) {
			int sign = column->isInCurrentSolution() ? -1 : 1;
			double cost_ = column->getCost();

			IloNumColumn col = cost(cost_ * sign);
			for (auto contrib : column->getContribs()) {
				if (contrib.second != 0) {
					col += constraints[constraintsIds[contrib.first]](sign* contrib.second);
				}
			}

			col += constraints[n_tasks](normalizationConstraint_[i]);

			/*
			if (acSupport_.count(i)) {
				col += constraints[constraintsIndexesAC[i]](-1);
			}
			else if (i == acolId_) {
				for (int j = 0; j < acSupport_.size(); j++) {
					col += constraints[n_tasks + 1 + j](1);
				}
			} */
			
			colsIndices.push_back(i);
			vars.add(IloNumVar(col, 0, IloInfinity));

			varsIds[i] = vars.getSize() - 1;


			col.end();
		}
	}

	cplex.setParam(IloCplex::Param::RootAlgorithm, 2);
	cplex.setParam(IloCplex::Param::Threads, 8);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Param::Advance, 1);
}

double ComplementaryProblemDis::solve(std::vector<double>* solution, std::map<int, int> originalProblemColumns,
std::map<std::string, std::pair<int, int>> tasksMapping, std::vector<double>& duals, bool first, bool export_) {

	std::cout << vars.getSize() << " variables dans le probleme complementaire desagrege." << std::endl;
	// R�solution du probl�me
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	int n = 1;
	if (originalProblemColumns.size()) {
		//cplex.setParam(IloCplex::Param::Simplex::Limits::Iterations, n);
	}
	
	IloNumArray duals_cplex(env, constraints.getSize());
	for(int i = 0; i < constraints.getSize()-1;i++) {
		duals_cplex[i] = duals[tasksMapping[psolutionMethod_->tasks_[i]].first];
	}

	duals_cplex[constraints.getSize()-1] = duals[duals.size() - 1];
	IloNumArray vals2;
	IloNumVarArray vars2;
	if(first) {
		cplex.setStart(vals2, 0, vars2, 0, duals_cplex, constraints);
	}


	bool success = cplex.solve();
	if (originalProblemColumns.size()) {
		while (cplex.getStatus() != IloAlgorithm::Status::Optimal) {
			if (true) {

				std::set<int> oc;
				std::set<int> goodvars;
				for (int i = 0; i < vars.getSize(); i++) {
					if (cplex.getBasisStatus(vars[i]) == IloCplex::BasisStatus::Basic) {
						if (!oc.count(originalProblemColumns[colsIndices[i]])) {
							goodvars.insert(i);
						}
						oc.insert(originalProblemColumns[colsIndices[i]]);
						
					}
				}


				for (int i = 0; i < vars.getSize(); i++) {
					if (!oc.count(originalProblemColumns[colsIndices[i]]) || goodvars.count(i)) {
						vars[i].setBounds(0, IloInfinity);
					}
					else {
						vars[i].setBounds(0, 0);
					}
				}

				//cplex.setParam(IloCplex::Param::Simplex::Limits::LowerObj, objective-1e-3);
			}
			//n += 1;
			cplex.setParam(IloCplex::Param::Simplex::Limits::Iterations, n);
			success = cplex.solve();

		}

		std::cout << "Boucle finie " << std::endl;
	}
	double objective;
	std::set<int> coveredTasks;
	if (success) {
		
		IloNumArray vals(env);
		cplex.getValues(vars, vals);
		double mean = 0;
		double newObjective = 0;
		int n = 0;
		for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
			solution->push_back(0);
		}

		for (int i = 0; i < vals.getSize(); i++) {
			if (vals[i] > 1e-5 && vars[i].getUB() > 0) {
				std::cout << originalProblemColumns[colsIndices[i]] << std::endl;
				mean += vals[i];
				n += 1;
				int sign = psolutionMethod_->columns_[colsIndices[i]]->isInCurrentSolution() ? -1 : 1;
				newObjective += sign * psolutionMethod_->columns_[colsIndices[i]]->getCost() * vals[i];
				for (auto contrib : psolutionMethod_->columns_[colsIndices[i]]->getContribs()) {
					coveredTasks.insert(tasksMapping[contrib.first].first);
				}
			}
			solution->at(colsIndices[i]) = vals[i];
		}

		objective = newObjective;
		int n2 = 0;
		std::set<std::string> constraintsToRemove;
		// Suppression des colonnes
		for (int i = 0; i < vars.getSize(); i++) {
			bool delete_ = false;
			for (auto contrib : psolutionMethod_->columns_[colsIndices[i]]->getContribs()) {
				if (contrib.second != 0 && coveredTasks.count(tasksMapping[contrib.first].first)) {
					vars[i].setBounds(0, 0);
					delete_ = true;
					if(!constraintsToRemove.count(contrib.first)) {
						constraintsToRemove.insert(contrib.first);
					}
				}
			}

			if (delete_) {
				n2 += 1;
			}
		}

		std::cout << n2 << " variables mises à 0" << std::endl;
		// Suppression des contraintes
		/*for (auto task : constraintsToRemove) {
			mod.remove(constraints[constraintsIds[task]]);
		}*/

	}
	else {
		return 10;
	}

	return objective;
}

void ComplementaryProblemDis::destroy() {
	env.end();
}
