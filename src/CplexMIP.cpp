#include "CplexMIP.h"

// Solve CPLEX MIP
double CplexMIP::solve(std::vector<int>& currentSolution, std::vector<int>* solution, bool relaxation, std::string path, std::string export_path) {
	//Construction du probl�me
	IloEnv env;
	IloModel mod(env);
	IloNumVarArray vars(env);

	//Valeurs initiales
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	std::ofstream cplexOut;

	std::vector<std::string> activeConstraints_ = psolutionMethod_->tasks_;
	std::map<std::string, int> tasksIndexes;
	for (int i = 0; i < activeConstraints_.size(); i++) {
		tasksIndexes[psolutionMethod_->tasks_[i]] = i;
	}
	// D�finition des contraintes
	IloRangeArray constraints(env, activeConstraints_.size(), 0, 0);
	std::cout << constraints.getSize() << "|" << activeConstraints_.size() << std::endl;
	for (int i = 0; i < activeConstraints_.size(); i++) {
		constraints[i].setBounds(psolutionMethod_->rhs_[i], psolutionMethod_->rhs_[i]);
		mod.add(constraints[i]);
	}

	IloObjective cost = IloAdd(mod, IloMinimize(env));
	cost.setConstant(psolutionMethod_->fixed_cost_);
	std::cout << cost.getConstant() << " is the fixed cost." << std::endl;

	std::vector<int> varsIds;

	// Rajout des variables

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {

		IB_Column* column = psolutionMethod_->columns_[i];
		if (true) {
			IloNumColumn col = cost(column->getCost());
			for (auto contrib : column->getContribs()) {
				col += constraints[tasksIndexes[contrib.first]](contrib.second);
			}

			if (!relaxation) {
				vars.add(IloBoolVar(col));
			}
			else {
				vars.add(IloNumVar(col, 0, 1));
			}
			if (currentSolution.size() > 0 && !relaxation) {
				startVar.add(vars[vars.getSize() - 1]);
				startVal.add((currentSolution)[i]);
			}
			
			varsIds.push_back(i);

			col.end();
		}
	}

	std::cout << vars.getSize() << " variables" << std::endl;
	// R�solution du probl�me
	IloCplex cplex(mod);
	if (export_path != "") {
		std::cout << "oui" << std::endl;
		cplex.exportModel(export_path.c_str());
		env.end();
		return 0;
	}
	cplex.setParam(IloCplex::Param::Threads, 8);
	if (currentSolution.size() > 0 && !relaxation) {
		cplex.addMIPStart(startVar, startVal);
	}
	if (path != "") {
		cplexOut.open(path + "/cplex.txt");
		cplex.setOut(cplexOut);
	}
	//cplex.setParam(IloCplex::PreInd, 0);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.005);
	cplex.setParam(IloCplex::Param::TimeLimit, 3600 * 5);

	bool success = cplex.solve();
	double objective;
	if (success) {
		IloNumArray vals(env);
		cplex.getValues(vars, vals);
		for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
			solution->push_back(0);
		}

		for (int i = 0; i < vals.getSize(); i++) {
			solution->at(varsIds[i]) = vals[i] > 1e-4 ? 1:0;
		}

		objective = cplex.getObjValue();
		double bestBound = objective * (1-cplex.getMIPRelativeGap());

		if (path != "") {
			std::ofstream cplexSolution;
			cplexSolution.open(path + "/cplex_solution.txt");
			cplexSolution << objective << " " << bestBound;
			cplexSolution.close();
		}

	}
	else {
		cplexOut.close();
		throw std::exception();
	}

	env.end();

	if (path != "") {
		cplexOut.close();
	}
	return objective;
}

// Solve CPLEX MIP from path
double CplexMIP::solveFromFile(std::string path, std::string solution_file) {
	IloEnv env;
	IloModel mod(env);
	IloCplex cplex(env);
	IloObjective obj;
	IloRangeArray constraints(env);
	IloNumVarArray vars(env);

	cplex.importModel(mod, path.c_str(), obj, vars, constraints);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.005);
	cplex.setParam(IloCplex::Param::TimeLimit, 3600 * 5);
	cplex.extract(mod);
	cplex.solve();
	IloNumArray vals(env);
	cplex.getValues(vars, vals);

	std::ofstream solution;
	solution.open(solution_file);
	for (int i = 0; i < vars.getSize(); i++) {
		solution << vars[i].getName() << " " << vals[i] << std::endl;
	}
	solution.close();

	return cplex.getObjValue();
	
}
