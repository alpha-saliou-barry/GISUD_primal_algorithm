#include "ISUD.h"
#include "IB_ComplementaryProblem.h"
#include "IB_CompatibilityChecker.h"
#include "IB_ReducedCP.h"
#include "IB_ReducedProblem.h"
#include "CplexMIP.h"
#include <cmath>
#include <iomanip>
#include <chrono>
#include <thread>
#include "ComplementaryProblemDis.h"
#include "DisComputer.h"

// Compute incompatibility degrees of the columns 'columns'
void calcIncompatibilityDegrees(IncompatibilityDegree id, std::vector<IB_Column *> columns)
{

	for (auto column : columns)
	{
		std::set<std::string> involvedColumns;
		column->setPhase(id.getIncompatibilityDegree(column, NULL, &involvedColumns));
		column->setCompatibleby(involvedColumns);
	}
}

// Constructor of ISUD
ISUD::ISUD(ISUD_Base *problem, bool addColumns, bool checkBinaryCompatibility, bool disE, bool compete_) : psolutionMethod_(problem), bcompatibilityChecker_(problem), addColumns_(addColumns), checkBinaryCompatibility_(checkBinaryCompatibility), compete(compete_)
{
	disEnabled = disE;
	currentCost_ = problem->fixed_cost_;
	for (int i = 0; i < problem->columns_.size(); i++)
	{
		IB_Column *column = problem->columns_[i];
		if (column->isInCurrentSolution())
		{
			currentSolution_.push_back(1);
			currentCost_ += column->getCost();
		}
		else
		{
			currentSolution_.push_back(0);
		}
	}

	gap = 0.05;
}

// Search for an integer direction in the support of one direction
bool ISUD::searchSubDirection(std::vector<int> in_columns, std::vector<int> *out_columns, std::vector<int> *colsIn)
{
	IloEnv env;
	IloModel mod(env);
	IloNumVarArray vars(env);

	int n_tasks = psolutionMethod_->tasks_.size();

	// D�finition des contraintes
	IloRangeArray constraints(env, n_tasks + 1, 0, 0);
	constraints[constraints.getSize() - 1].setBounds(1, IloInfinity);
	for (int i = 0; i < constraints.getSize(); i++)
	{
		mod.add(constraints[i]);
	}

	IloObjective cost = IloAdd(mod, IloMinimize(env));
	std::vector<int> colsIndices;

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];

		if (column->isInP())
		{
			int sign = column->isInCurrentSolution() ? -1 : 1;
			IloNumColumn col = cost(column->getCost() * sign);
			for (int j = 0; j < n_tasks; j++)
			{
				int contrib = column->findContribution(psolutionMethod_->tasks_[j]);
				if (contrib != 0)
				{
					col += constraints[j](sign * contrib);
				}
			}
			col += constraints[constraints.getSize() - 1](1);
			colsIndices.push_back(i);
			vars.add(IloBoolVar(col));
			col.end();
		}
	}

	for (int i = 0; i < in_columns.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[in_columns[i]];
		std::cout << "in_column  : " << std::endl;
		int sign = column->isInCurrentSolution() ? -1 : 1;
		IloNumColumn col = cost(column->getCost() * sign);
		for (int j = 0; j < n_tasks; j++)
		{
			int contrib = column->findContribution(psolutionMethod_->tasks_[j]);
			if (contrib != 0)
			{
				col += constraints[j](sign * contrib);
			}
		}

		col += constraints[constraints.getSize() - 1](1);
		colsIndices.push_back(in_columns[i]);
		vars.add(IloBoolVar(col));

		col.end();
	}

	std::cout << vars.getSize() << " variables dans le probleme de sous direction." << std::endl;
	// R�solution du probl�me
	IloCplex cplex(mod);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::Simplex::Display, 0);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setOut(env.getNullStream());

	bool success = cplex.solve();
	std::cout << (success ? cplex.getObjValue() : 0) << std::endl;
	if (!success || (cplex.getObjValue() >= -1e-4))
	{
		return false;
	}

	IloNumArray vals(env);
	cplex.getValues(vars, vals);

	for (int i = 0; i < vals.getSize(); i++)
	{
		if (vals[i] > 1e-5 && i < colsIndices.size())
		{
			if (psolutionMethod_->columns_[colsIndices[i]]->isInCurrentSolution())
			{
				out_columns->push_back(colsIndices[i]);
			}
			else
			{
				colsIn->push_back(colsIndices[i]);
			}
		}
	}

	return true;
}

// Add a line to the output
void ISUD::addLine(double newCost, double amelioration, int n_iterations, int n_added_columns, int n_success_add_columns, bool has_zoom, int phaseMax, bool addedColumn, int n_added_columns_it, int pivot_distance, double iteration_time, double global_time)
{
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << n_iterations;
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << pivot_distance;
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << std::fixed << std::setprecision(6) << (global_time / 1000);
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << std::fixed << std::setprecision(6) << (iteration_time / 1000);
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << std::fixed << std::setprecision(8) << newCost;
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << std::fixed << std::setprecision(8) << amelioration;
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << phaseMax;
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << (has_zoom ? "Oui" : "Non");
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << (addedColumn ? "Oui" : "Non");
	if (addedColumn && !disEnabled)
	{
		file << " (" << n_added_columns_it << ")";
	}
	if (n_added_columns > 0)
	{
		file << std::left << std::setw(nameWidth) << std::setfill(separator) << std::fixed << std::setprecision(4) << ((((double)n_success_add_columns) / ((double)n_added_columns)) * 100.0);
	}
	else
	{
		file << std::left << std::setw(nameWidth) << std::setfill(separator) << "_";
	}

	if (addedColumn && disEnabled)
	{
		file << std::left << std::setw(nameWidth) << dis_problem_size;
	}

	file << std::endl;
}

//Return if the support of one direction can be included in an integer direction
bool ISUD::canBeInIntegerDirection(IB_Column *column, std::vector<double> &solution, std::vector<int> *support, int acolId)
{
	IB_CompatibilityChecker cc(psolutionMethod_);
	Eigen::VectorXf newColumn(psolutionMethod_->tasks_.size());
	for (int i = 0; i < newColumn.size(); i++)
	{
		newColumn(i) = 0;
	}

	Eigen::VectorXf b(psolutionMethod_->tasks_.size());
	for (int i = 0; i < psolutionMethod_->rhs_.size(); i++)
	{
		b(i) = psolutionMethod_->rhs_[i];
	}

	double cost = 0;

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (i != acolId)
		{
			if (solution[i] > 1e-5 || (support != NULL && std::find(support->begin(), support->end(), i) != support->end()))
			{
				int sign = column->isInCurrentSolution() ? 1 : -1;
				newColumn += sign * cc.columnToEigenVector(column);
				if (sign == -1)
				{
					b -= cc.columnToEigenVector(column);
				}

				cost += sign * column->getCost();
			}
		}
	}

	if ((b.array() > -1e-4).all())
	{
		std::vector<std::pair<std::string, int>> contribs;
		for (int i = 0; i < newColumn.size(); i++)
		{
			if (fabs(newColumn(i)) > 1e-4)
			{
				contribs.push_back(std::pair<std::string, int>(psolutionMethod_->tasks_[i], (int)round(newColumn(i))));
			}
		}

		column->setContribs(contribs);
		column->reorderContribs(psolutionMethod_->tasks_);
		column->setCost(cost);
		return true;
	}

	return false;
}

// Return binary compatible column of negative reduced cost
bool ISUD::getBCompatibleColumn(std::vector<int> *colsIn, std::vector<int> *colsOut)
{
	int bcColumn = -1;
	double bcCost = 0;

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (!column->isInCurrentSolution() && column->getState() == IB_Column::STATE_COMPATIBLE)
		{
			double reducedCost = column->getCost();

			for (int colId : bcompatibilityChecker_.makeCompatibles_[i])
			{
				reducedCost -= psolutionMethod_->columns_[colId]->getCost();
			}

			if (reducedCost < -1e-4)
			{
				if (bcColumn == -1 || reducedCost < bcCost)
				{
					bcColumn = i;
					bcCost = reducedCost;
				}
			}
		}
	}

	if (bcColumn == -1)
	{
		return false;
	}

	int i = bcColumn;
	colsIn->push_back(i);
	for (int colId : bcompatibilityChecker_.makeCompatibles_[i])
	{
		colsOut->push_back(colId);
	}

	return true;
}

// Return binary compatible column of negative reduced cost with incompatibility degree
bool ISUD::getBCompatibleColumnId(std::vector<int> *colsIn, std::vector<int> *colsOut)
{
	int bcColumn = -1;
	double bcCost = 0;

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (!column->isInCurrentSolution() && column->getPhase() == 0)
		{
			double reducedCost = column->getCost();

			for (std::string colName : column->getCompatibleBy())
			{
				int colId = psolutionMethod_->getColumnIndice(colName);
				reducedCost -= psolutionMethod_->columns_[colId]->getCost();
			}

			if (reducedCost < -1e-4)
			{
				if (bcColumn == -1 || reducedCost < bcCost)
				{
					bcColumn = i;
					bcCost = reducedCost;
				}
			}
		}
	}

	if (bcColumn == -1)
	{
		return false;
	}

	{
		int i = bcColumn;
		colsIn->push_back(i);
		for (std::string colName : psolutionMethod_->columns_[bcColumn]->getCompatibleBy())
		{
			int colId = psolutionMethod_->getColumnIndice(colName);
			colsOut->push_back(colId);
		}
	}

	return true;
}

// Pivot columns colsIn, colsOut in the solution, recompute compatibilities if recomputeCompatibilities is true
void ISUD::pivotColumnsInSolution(std::vector<int> &colsIn, std::vector<int> &colsOut, bool recomputeCompatibilities)
{
	std::set<std::string> coveredTasks;

	for (int colId : colsIn)
	{
		std::cout << "ColIn : " << colId << std::endl;
		IB_Column *column = psolutionMethod_->columns_[colId];
		column->setInCurrentSolution();
		currentSolution_[colId] = 1;
		currentCost_ += column->getCost();
		for (auto contrib : column->getContribs())
		{
			if (contrib.second != 0)
			{
				coveredTasks.insert(contrib.first);
			}
		}
	}

	for (int colId : colsOut)
	{
		std::cout << "ColOut : " << colId << std::endl;
		IB_Column *column = psolutionMethod_->columns_[colId];
		column->setOutCurrentSolution();
		column->setIncompatible();
		currentSolution_[colId] = 0;
		currentCost_ -= column->getCost();

		for (auto contrib : column->getContribs())
		{
			if (contrib.second != 0)
			{
				coveredTasks.insert(contrib.first);
			}
		}
	}

	double cost = psolutionMethod_->fixed_cost_;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		if (psolutionMethod_->columns_[i]->isInCurrentSolution())
		{
			cost += psolutionMethod_->columns_[i]->getCost();
		}
	}

	currentCost_ = cost;

	std::cout << "Nouveau cout : " << currentCost_ << std::endl;

	// On recalcule les compatibilit�s binaires des colonnes
	if (checkBinaryCompatibility_)
	{
		std::vector<int> columnsToRecomputeCompatibility;
		for (int colId : colsOut)
		{
			columnsToRecomputeCompatibility.push_back(colId);
		}

		std::set<int> columnsOut;
		for (int colId : colsOut)
		{
			columnsOut.insert(colId);
		}

		std::set<int> columnsIn;
		for (int colId : colsIn)
		{
			columnsIn.insert(colId);
		}

		for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
		{
			IB_Column *column = psolutionMethod_->columns_[i];

			bool recompute = false;
			if (!column->isInCurrentSolution() && column->getState() == IB_Column::STATE_COMPATIBLE)
			{
				for (auto colId : bcompatibilityChecker_.makeCompatibles_[i])
				{
					if (columnsOut.count(colId))
					{
						recompute = true;
						break;
					}
				}
			}
			else if (!column->isInCurrentSolution() && column->getState() == IB_Column::STATE_INCOMPATIBLE)
			{
				for (auto colId : bcompatibilityChecker_.subsetColumns_[i])
				{
					if (columnsIn.count(colId))
					{
						recompute = true;
					}
				}
			}

			if (recompute)
			{
				columnsToRecomputeCompatibility.push_back(i);
			}
		}

		std::cout << columnsToRecomputeCompatibility.size() << " colonnes dont on recalcule la compatibibilite binaire" << std::endl;
		for (auto colId : columnsToRecomputeCompatibility)
		{
			1 + 1;
			bcompatibilityChecker_.updateCompatibilityStatus(colId);
		}
	}

	std::cout << "Recalcul des degr�s d'incompatibilit�s" << std::endl;
	std::vector<IB_Column *> positiveColumns;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		if (psolutionMethod_->columns_[i]->isInCurrentSolution())
		{
			positiveColumns.push_back(psolutionMethod_->columns_[i]);
		}
	}

	// On recalcule les r d'incompatibilit�s
	if (recomputeCompatibilities)
	{

		std::vector<IB_Column *> columns_to_recompute;
		std::set<std::string> colsDones;
		for (std::string taskName : coveredTasks)
		{
			if (psolutionMethod_->columnsPerTask.count(taskName))
			{
				for (auto column : psolutionMethod_->columnsPerTask[taskName])
				{
					if (!colsDones.count(column))
					{
						colsDones.insert(column);
						columns_to_recompute.push_back(psolutionMethod_->columns_[psolutionMethod_->colsIndices[column]]);
					}
				}
			}
		}

		int n_threads = 8;

		if (columns_to_recompute.size() <= 100)
		{
			n_threads = 1;
		}

		std::vector<std::vector<IB_Column *>> thread_tasks;
		std::vector<IB_Column *> current_vector;
		for (int i = 0; i < columns_to_recompute.size(); i++)
		{
			if (current_vector.size() >= ((float)columns_to_recompute.size()) / n_threads)
			{
				thread_tasks.push_back(current_vector);
				current_vector = {};
			}

			current_vector.push_back(columns_to_recompute[i]);
		}

		if (current_vector.size())
		{
			thread_tasks.push_back(current_vector);
		}

		std::cout << thread_tasks.size() << " est la taille." << std::endl;
		std::vector<IB_Column *> positiveColumns;
		for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
		{
			if (psolutionMethod_->columns_[i]->isInCurrentSolution())
			{
				positiveColumns.push_back(psolutionMethod_->columns_[i]);
			}
		}

		IncompatibilityDegree id(psolutionMethod_, positiveColumns, psolutionMethod_->tasks_);
		std::vector<std::thread> threads_;
		for (auto thread_task : thread_tasks)
		{
			threads_.push_back(std::thread(calcIncompatibilityDegrees, id, thread_task));
		}

		for (int i = 0; i < threads_.size(); i++)
		{
			threads_[i].join();
		}
	}
	else
	{
		for (std::string taskName : coveredTasks)
		{
			if (psolutionMethod_->columnsPerTask.count(taskName))
			{
				for (auto column : psolutionMethod_->columnsPerTask[taskName])
				{
					psolutionMethod_->columns_[psolutionMethod_->colsIndices[column]]->setPhase(-1);
				}
			}
		}
	}

	IB_CompatibilityChecker ib(psolutionMethod_);

	Eigen::VectorXf sum(psolutionMethod_->rhs_.size());
	for (int i = 0; i < psolutionMethod_->rhs_.size(); i++)
	{
		sum(i) = psolutionMethod_->rhs_[i];
	}

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (column->isInCurrentSolution())
		{
			sum -= ib.columnToEigenVector(column);
		}
	}

	if ((sum.cast<int>().array() != 0).any())
	{
		std::cout << sum.transpose() << std::endl;
		throw std::exception();
	}
}

// Return if the solution is integral
bool ISUD::isIntegral(std::vector<double> &solution)
{
	std::vector<double> nonZeroComponents;
	for (int i = 0; i < solution.size(); i++)
	{
		if (solution[i] > 1e-5)
		{
			nonZeroComponents.push_back(solution[i]);
		}
	}

	if (!nonZeroComponents.size())
	{
		return false;
	}

	double commonValue = nonZeroComponents[0];
	for (int i = 1; i < nonZeroComponents.size(); i++)
	{
		if (fabs(nonZeroComponents[i] - commonValue) > 1e-3)
		{
			return false;
		}
	}

	return true;
}

// ZOOM procedure
bool ISUD::zoom(int isudPhase, std::vector<int> seqPhases, std::vector<double> &solution, std::vector<int> *colsIn, std::vector<int> *colsOut)
{
	seqPhases = {isudPhase};
	std::vector<int> rpPhases = {isudPhase};
	std::cout << "Procedure zoom" << std::endl;
	bool cpSuccess = true;
	int n_iterations = 0;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (column->isInCurrentSolution() || solution[i] > 1e-5)
		{
			column->changeInP(true);
		}
		else
		{
			column->changeInP(false);
		}
	}
	while (cpSuccess)
	{
		std::cout << solution.size() << " taille de la solution." << std::endl;
		std::vector<double> incompatibleCosts;
		std::vector<int> incompatibleColsIndices;
		std::vector<Eigen::VectorXf> incompatibleVecs;
		std::vector<bool> compatibilityCheck;
		std::vector<std::string> activeConstraintsCP, activeConstraintsRP;
		std::set<int> activePColumns;
		std::vector<double> dualVariables;
		if (true)
		{
			std::vector<IB_Column *> positiveColumns;
                        std::vector<IB_Column*> zeroColumns;
			std::vector<int> positiveColumnsIndices;
			for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
			{
				if (psolutionMethod_->columns_[i]->isInP())
				{
					positiveColumns.push_back(psolutionMethod_->columns_[i]);
					positiveColumnsIndices.push_back(i);
				} else {
                                  if(isudPhase == -1 || psolutionMethod_->columns_[i]->getPhase() <= isudPhase) {
                                       zeroColumns.push_back(psolutionMethod_->columns_[i]);
                                  } else {
                                       psolutionMethod_->columns_[i]->setIncompatible();
                                  }
                                }
			}

			// Degr� d'incompatibilit�

			// Construction du RP et du CP
			IB_CompatibilityChecker compatibilityChecker(psolutionMethod_);
			compatibilityChecker.calcIndependentMatrix();
			compatibilityChecker.calcInverseStructure();

			std::vector<IB_Column *> incompatibleColumns;

			int rp_normal_size = 0;
			compatibilityChecker.isLinearlyCompatible(zeroColumns, compatibilityCheck);
			// std::cout << "Check des compatibilit�s binaires" << std::endl;
			for (int i = 0; i < zeroColumns.size(); i++)
			{
				IB_Column *column = zeroColumns[i];
				if (!column->isInP())
				{
					std::vector<int> columns_compatible;
					if (compatibilityCheck[i])
					{
						column->setCompatibleLinearly();
						rp_normal_size += 1;
					}
					else
					{
						column->setIncompatible();
						incompatibleColumns.push_back(column);
						incompatibleColsIndices.push_back(i);
					}
				}
				else
				{
					rp_normal_size += 1;
				}
			}

			// R�solution du probl�me r�duit

			bool solveAgainrp = true;
			bool changeConstraints = false;
			std::vector<std::string> activeConstraints2 = compatibilityChecker.getActiveConstraints();
			double gapValue = 0.005;
			int rp_max_size = psolutionMethod_->columns_.size() * 100 / 100;
			int rpPhase = 0;
			while (solveAgainrp)
			{
				std::vector<std::string> ac = compatibilityChecker.getActiveConstraints();
				// std::vector<std::string> ac = psolutionMethod_->tasks_;
				IB_ReducedProblem rp(psolutionMethod_, ac, gapValue);
				std::vector<int> newSolution;

				solveAgainrp = false;
				double objective = rp.solveProblem(currentSolution_, &newSolution, currentCost_, rp_max_size);
				std::cout << currentCost_ - objective - 1e-4 << std::endl;
				if (objective < currentCost_ - 1e-4)
				{
					Eigen::VectorXf sum(psolutionMethod_->rhs_.size());
					for (int i = 0; i < psolutionMethod_->rhs_.size(); i++)
					{
						sum(i) = psolutionMethod_->rhs_[i];
					}

					for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
					{
						IB_Column *column = psolutionMethod_->columns_[i];
						if (newSolution[i] > 1e-4)
						{
							sum -= compatibilityChecker.columnToEigenVector(column);
						}
					}

					if ((sum.cast<int>().array() != 0).any())
					{
						solveAgainrp = true;
						std::cout << "SOLUTION DU RP invalide, on recommence" << std::endl;
						changeConstraints = true;
						for (int i = 0; i < sum.size(); i++)
						{
							if (fabs(sum(i)) > 1e-4)
							{
								compatibilityChecker.addActiveConstraint(i);
							}
						}
					}
					else
					{
						std::cout << "Nouvelle solution avec le RP" << std::endl;
						for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
						{
							IB_Column *column = psolutionMethod_->columns_[i];

							if (column->isInCurrentSolution() && newSolution[i] < 1e-4)
							{
								colsOut->push_back(i);
							}
							else if (!column->isInCurrentSolution() && newSolution[i] > 1 - 1e-4)
							{
								colsIn->push_back(i);
							}
						}

						return true;
					}
				}
				else
				{
					if (rp_normal_size <= rp_max_size)
					{
						rpPhase = 0;
						if (gapValue > 0.005 + 1e-4)
						{
							rp_max_size = psolutionMethod_->columns_.size() * 100 / 100;
							gapValue = 0.005;
							solveAgainrp = true;
						}
					}
					else
					{
						rp_max_size += psolutionMethod_->columns_.size() * 3 / 100;
						solveAgainrp = true;
					}
					std::cout << "Constance de la solution" << std::endl;

					if (!solveAgainrp)
					{
						rp.getDuals(&dualVariables);
						activeConstraintsRP = ac;
					}
				}
			}

			if (changeConstraints)
			{
				compatibilityChecker.setActiveConstraints(activeConstraints2);
			}
			// wcompatibilityChecker.getIncompatibleVectorAndCost(incompatibleColumns, &incompatibleVecs, &incompatibleCosts);
			activeConstraintsCP = compatibilityChecker.getCPActiveConstraints();
			activePColumns = compatibilityChecker.getActiveColumns();
		}

		std::vector<Eigen::VectorXf *> incompatibleVecs2;
		for (int i = 0; i < incompatibleVecs.size(); i++)
		{
			incompatibleVecs2.push_back(&incompatibleVecs[i]);
		}

		for (int i = 0; i < seqPhases.size(); i++)
		{
			int phase = seqPhases[i];
			std::vector<double> incompatibleCosts_phase;
			std::vector<int> incompatibleColsIndices_phase;
			std::vector<Eigen::VectorXf *> incompatibleVecs2_phase;
			for (int j = 0; j < incompatibleVecs2.size(); j++)
			{
				if (phase == -1 || psolutionMethod_->columns_[incompatibleColsIndices[j]]->getPhase() <= phase)
				{
					incompatibleCosts_phase.push_back(incompatibleCosts[j]);
					incompatibleColsIndices_phase.push_back(incompatibleColsIndices[j]);
					incompatibleVecs2_phase.push_back(incompatibleVecs2[j]);
				}
			}

			IB_ReducedCP reducedCP(psolutionMethod_, activePColumns, activeConstraintsCP, incompatibleColsIndices_phase, incompatibleVecs2_phase, incompatibleCosts_phase);
			std::vector<IB_Column *> colsToAdd;
			std::vector<int> colsToAddIndices;
			std::cout << "Resolution du probleme complementaire reduit en phase " << phase << std::endl;

			cpSuccess = reducedCP.solve(currentCost_, dualVariables, activeConstraintsRP, &colsToAdd, &colsToAddIndices, 3, phase);
			if (cpSuccess && colsToAdd.size() > 0)
			{
				int n = 0;
				// On v�rifie si la solution est enti�re :
				std::set<std::string> coveredTasks;
				bool column_disjoint = true;
				for (auto col : colsToAdd)
				{
					for (auto task_pair : col->getContribs())
					{
						if (coveredTasks.count(task_pair.first))
						{
							column_disjoint = false;
						}

						coveredTasks.insert(task_pair.first);
					}
				}
				if (column_disjoint)
				{
					// On cherche une sous direction de co�t r�duit n�gatif
					if (searchSubDirection(colsToAddIndices, colsOut, colsIn))
					{
						return true;
					}
				}
				for (auto col : colsToAdd)
				{
					col->changeInP(true);
					n += 1;
				}
				std::cout << n << " colonnes rajoutees dans p" << std::endl;
				break;
			}
		}

		if (!cpSuccess)
		{
			std::cout << "Echec du CP, fin de la resolution." << std::endl;
			return false;
		}

		n_iterations += 1;
	}

	if (cpSuccess)
	{
		std::cout << "Echec de la procedure zoom, fin." << std::endl;
		return false;
	}

	return false;
}

// Complementary problem with column addition
std::pair<bool, int> ISUD::cpWithArtificialColumn(int acolId, std::set<int> colsIn, std::set<int> colsOut,
												  std::vector<int> &ncolsIn, std::vector<int> &ncolsOut, int initial_phase, int n_a_cols, double previous_objective)
{
	std::vector<int> support;
	std::set<std::string> coveredTasks;
	for (auto colId : colsOut)
	{
		for (auto contrib : psolutionMethod_->columns_[colId]->getContribs())
		{
			if (contrib.second != 0)
			{
				coveredTasks.insert(contrib.first);
			}
		}
	}

	std::cout << "Rajout d'une colonne : " << psolutionMethod_->columns_[acolId]->getName() << std::endl;
	psolutionMethod_->columns_[acolId]->setInCurrentSolution();

	for (auto colId : colsIn)
	{
		// psolutionMethod_->columns_[colId]->setInCurrentSolution();
		support.push_back(colId);
	}

	for (auto colId : colsOut)
	{
		// psolutionMethod_->columns_[colId]->setOutCurrentSolution();
		support.push_back(colId);
	}

	std::vector<IB_Column *> positiveColumns;
	std::vector<int> pastIds;

	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (column->isInCurrentSolution())
		{
			positiveColumns.push_back(column);
		}

		pastIds.push_back(column->getPhase());
	}

	/*
	std::cout << "Calcul des degres d'incompatibilites" << std::endl;

	IncompatibilityDegree id(psolutionMethod_, positiveColumns, psolutionMethod_->tasks_);
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		IB_Column* column = psolutionMethod_->columns_[i];
		if (!column->isInCurrentSolution()) {
			bool recompute = false;
			for (auto contribPair : column->getContribs()) {
				if (contribPair.second != 0 && coveredTasks.count(contribPair.first)) {
					recompute = true;
					break;
				}
			}

			if (recompute) {
				column->setPhase(id.getIncompatibilityDegree(column, NULL));
			}
		}
	}*/

	std::vector<int> phaseSeq;
	for (auto phase : getPhaseSequence())
	{
		if (phase >= initial_phase || phase == -1)
		{
			phaseSeq.push_back(phase);
		}
	}

	std::vector<double> solution;

	for (int i = 0; i < phaseSeq.size(); i++)
	{
		std::cout << "Resolution du probleme complementaire en phase : " << phaseSeq[i] << std::endl;
		solution.clear();
		std::vector<int> support2;
		for (auto column_id : support)
		{
			IB_Column *column = psolutionMethod_->columns_[column_id];
			if (column->isInCurrentSolution() || (phaseSeq[i] == -1 || column->getPhase() <= phaseSeq[i]))
			{
				support2.push_back(column_id);
			}
		}
		IB_ComplementaryProblem cp(psolutionMethod_, phaseSeq[i], acolId, &support2, support.size());
		cp.setPhase(phaseSeq[i]);
		double objective = -1;
		std::cout << "Penalisation : " << previous_objective * (support.size() + 200) << std::endl;
		cp.constructProblem(true, previous_objective * (support.size() + 200));
		objective = cp.solve(&solution);
		cp.destroy();

		bool good_objective = objective < -1e-3;

		if (!good_objective)
		{
			std::cout << "Echec du probleme complementaire" << std::endl;
			for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
			{
				IB_Column *column = psolutionMethod_->columns_[j];
				if (colsIn.count(j))
				{
					column->setOutCurrentSolution();
				}
				else if (colsOut.count(j))
				{
					column->setInCurrentSolution();
				}

				column->setPhase(pastIds[j]);
			}

			psolutionMethod_->removeColumn(acolId);
			return std::pair<bool, int>(false, n_a_cols);
		}

		/*if (solution[acolId] < 1e-4) {
			std::cout << "Colonne artificielle pas dans la solution." << std::endl;
			continue;
		}*/

		if (good_objective)
		{

			if (solution[acolId] > 1e-5)
			{
				std::cout << "colonne artificielle dans la solution" << std::endl;
			}

			std::cout << objective << std::endl;
			artificial_column_id += 1;
			IB_Column *artificial_column = new IB_Column("artificial_column_" + std::to_string(artificial_column_id),
														 std::vector<std::pair<std::string, int>>());

			if (isIntegral(solution) && solution[acolId] > 1e-5)
			{
				std::cout << "Direction entiere trouvee " << std::endl;
				ncolsIn.clear();
				ncolsOut.clear();

				for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
				{
					IB_Column *column = psolutionMethod_->columns_[j];
					if (solution[j] > 1e-5 && j != acolId)
					{
						if (column->isInCurrentSolution())
						{
							if (!colsIn.count(j))
							{
								ncolsOut.push_back(j);
							}
						}
						else
						{
							if (!colsOut.count(j))
							{
								ncolsIn.push_back(j);
							}
						}
					}
				}

				if (solution[acolId] > 1e-5)
				{
					std::cout << "colonne artificielle dans la solution" << std::endl;
					for (auto colId : colsIn)
					{
						ncolsIn.push_back(colId);
						// ncolsIn.push_back(colId);
					}

					for (auto colId : colsOut)
					{
						ncolsOut.push_back(colId);
						// ncolsOut.push_back(colId);
					}
				}

				for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
				{
					IB_Column *column = psolutionMethod_->columns_[j];
					if (colsIn.count(j))
					{
						column->setOutCurrentSolution();
					}
					else if (colsOut.count(j))
					{
						column->setInCurrentSolution();
					}

					column->setPhase(pastIds[j]);
				}
				std::cout << "ok1" << std::endl;
				psolutionMethod_->removeColumn(acolId);
				std::cout << "ok2" << std::endl;

				return std::pair<bool, int>(true, n_a_cols);
			}
			else if (canBeInIntegerDirection(artificial_column, solution, &support, acolId) && solution[acolId] > 1e-5)
			{

				std::set<int> newColsIn, newColsOut;
				for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
				{
					if (solution[j] > 1e-5)
					{
						if (j != acolId)
						{
							if (psolutionMethod_->columns_[j]->isInCurrentSolution())
							{
								newColsOut.insert(j);
							}
							else
							{
								newColsIn.insert(j);
							}
						}
					}
				}

				for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
				{
					IB_Column *column = psolutionMethod_->columns_[j];
					if (colsIn.count(j))
					{
						column->setOutCurrentSolution();
					}
					else if (colsOut.count(j))
					{
						column->setInCurrentSolution();
					}

					column->setPhase(pastIds[j]);
				}

				for (int j : support)
				{
					if (psolutionMethod_->columns_[j]->isInCurrentSolution())
					{
						newColsOut.insert(j);
					}
					else
					{
						newColsIn.insert(j);
					}
				}

				psolutionMethod_->removeColumn(acolId);
				int colId = psolutionMethod_->addArtificialColumn(artificial_column);
				return cpWithArtificialColumn(colId, newColsIn, newColsOut, ncolsIn, ncolsOut, initial_phase, n_a_cols + 1, previous_objective);
			}
			else
			{
				for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
				{
					IB_Column *column = psolutionMethod_->columns_[j];
					if (colsIn.count(j))
					{
						column->setOutCurrentSolution();
					}
					else if (colsOut.count(j))
					{
						column->setInCurrentSolution();
					}

					column->setPhase(pastIds[j]);
				}
				psolutionMethod_->removeColumn(acolId);
				return std::pair<bool, int>(false, n_a_cols);
			}
		}
	}

	for (int j = 0; j < psolutionMethod_->columns_.size(); j++)
	{
		IB_Column *column = psolutionMethod_->columns_[j];
		if (colsIn.count(j))
		{
			column->setOutCurrentSolution();
		}
		else if (colsOut.count(j))
		{
			column->setInCurrentSolution();
		}

		column->setPhase(pastIds[j]);
	}

	psolutionMethod_->removeColumn(acolId);

	return std::pair<bool, int>(false, n_a_cols);
	;
}

// Add row for the comparison of ZOOM and column addition
void ISUD::addCompeteRow(double amelioration_rc, int time_rc, double amelioration_zoom, int time_zoom, int disp_size, double last_objective, double remaining)
{
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << currentCost_;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << amelioration_rc;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << time_rc;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << amelioration_zoom;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << time_zoom;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << disp_size;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << last_objective;
	fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << remaining;
	fileC << std::endl;
}

// Main procedure of ISUD, solve the problem and stock the output to path
void ISUD::solve(std::string path)
{
	std::cout << psolutionMethod_->columns_.size() << " colonnes." << std::endl;
	std::cout << psolutionMethod_->tasks_.size() << " taches." << std::endl;
	std::string final_path = "";
	if (disEnabled)
	{
		if (compete)
		{
			final_path = path + "/sortie_isud_dis_compete.txt";
		}
		else
		{
			final_path = path + "/sortie_isud_dis.txt";
		}
	}
	else if (addColumns_)
	{
		if (compete)
		{
			final_path = path + "/sortie_isud_compete.txt";
		}
		else
		{
			final_path = path + "/sortie_isud.txt";
		}
	}
	else
	{
		final_path = path + "/sortie_zoom.txt";
	}

	file.open(path != "" ? final_path : "sortie.txt");
	if (compete)
	{
		fileC.open(path + (addColumns_ ? "/competition.txt" : "/competition_dis.txt"));
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Cout courant";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Amelioration RC";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Temps RC";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Amelioration ZOOM";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Temps ZOOM";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Size DISP";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Objective";
		fileC << std::left << std::setw(nameWidth) << std::setfill(separator) << "Remaining";
		fileC << std::endl;
	}
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Iteration";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Pivot Distance";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Temps total (s)";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Temps iteration (s)";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Meilleure solution entiere";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Amelioration";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Phase max";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Zoom";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Rajout de colonne";
	file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Pourcentage rajout de colonne";
	if (disEnabled)
	{
		file << std::left << std::setw(nameWidth) << std::setfill(separator) << "Taille PBS";
	}
	file << std::endl;

	IB_CompatibilityChecker ib(psolutionMethod_);

	Eigen::VectorXf sum(psolutionMethod_->rhs_.size());
	for (int i = 0; i < psolutionMethod_->rhs_.size(); i++)
	{
		sum(i) = 0;
	}

	int n_pos_cols = 0;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		IB_Column *column = psolutionMethod_->columns_[i];
		if (column->isInCurrentSolution())
		{
			sum += ib.columnToEigenVector(column);
			n_pos_cols += 1;
		}
	}

	std::cout << psolutionMethod_->columns_.size() << " colonnes." << std::endl;
	std::cout << n_pos_cols << " colonnes positives." << std::endl;
	std::cout << sum.transpose() << std::endl;

	// Initialisation de la compatibilit� binaire
	if (checkBinaryCompatibility_)
	{
		std::cout << "Calcul des compatibilites binaires des colonnes : " << std::endl;
		bcompatibilityChecker_.init();
		std::cout << "Compatibilites binaires calculees." << std::endl;
	}

	current_gap = 1;
	auto global_start = std::chrono::high_resolution_clock::now();
	std::cout << "Calcul de la valeur de la relaxation lineaire." << std::endl;
	std::vector<int> initialSolution = getCurrentSolution();
	CplexMIP cplex(psolutionMethod_);
	std::vector<int> solution;
	double objective = cplex.solve(initialSolution, &solution, true);
	bestBound = objective;
	double bound = bestBound;
	std::cout << "Relaxation linéaire : " << bound << std::endl;
	// Calcul des degr�s d'incompatibilit�s
	std::vector<IB_Column *> positiveColumns;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
	{
		if (psolutionMethod_->columns_[i]->isInCurrentSolution())
		{
			positiveColumns.push_back(psolutionMethod_->columns_[i]);
		}
	}

	std::cout << "Calcul des degres d'incompatibilites des colonnes" << std::endl;
	std::vector<IB_Column *> columns_to_recompute = psolutionMethod_->columns_;
	int n_threads = 8;

	if (columns_to_recompute.size() <= 100)
	{
		n_threads = 1;
	}

	std::vector<std::vector<IB_Column *>> thread_tasks;
	std::vector<IB_Column *> current_vector;
	for (int i = 0; i < columns_to_recompute.size(); i++)
	{
		if (current_vector.size() >= ((float)columns_to_recompute.size()) / n_threads)
		{
			thread_tasks.push_back(current_vector);
			current_vector = {};
		}

		current_vector.push_back(columns_to_recompute[i]);
	}

	if (current_vector.size())
	{
		thread_tasks.push_back(current_vector);
	}

	std::cout << thread_tasks.size() << " est la taille." << std::endl;

	std::vector<std::thread> threads_;
	for (auto thread_task : thread_tasks)
	{
		IncompatibilityDegree id(psolutionMethod_, positiveColumns, psolutionMethod_->tasks_);
		threads_.push_back(std::thread(calcIncompatibilityDegrees, id, thread_task));
	}

	for (int i = 0; i < threads_.size(); i++)
	{
		threads_[i].join();
	}

	std::cout << "Degres d'incompatibilites calcules" << std::endl;

	// On commence les it�rations
	bool solved = false;
	int n_iterations = 0;
	int n_added_columns = 0;
	int n_success_add_columns = 0;

	bool skipPhase = false;
	int previousPhaseMax = -1;
	while (!solved)
	{
		bool has_zoom = false;
		int phaseMax = -1;
		bool columnAdded = false;
		int pivot_distance = 1;
		int n_added_columns_it = 1;
		double last_objective = -10000;
		auto iteration_start = std::chrono::high_resolution_clock::now();
		std::vector<int> colsIn, colsOut;
		// Pivotage des colonnes compatibles binaires jusqu'� que ce ne soit plus possible
		if (checkBinaryCompatibility_)
		{
			bool failed = false;
			while (!failed)
			{
				failed = !getBCompatibleColumn(&colsIn, &colsOut);
				if (!failed)
				{
					std::cout << "Une colonne compatible binaire trouvee" << std::endl;
					pivotColumnsInSolution(colsIn, colsOut);

					colsIn.clear();
					colsOut.clear();
				}
			}
		}
		else
		{
			bool continue_ = true;
			while (continue_)
			{

				bool failed = false;
				bool has_compatible_column = false;
				while (!failed)
				{
					failed = !getBCompatibleColumnId(&colsIn, &colsOut);
					if (!failed)
					{
						has_compatible_column = true;
						std::cout << "Une colonne compatible binaire trouvee" << std::endl;
						pivotColumnsInSolution(colsIn, colsOut, false);

						colsIn.clear();
						colsOut.clear();
					}
				}

				positiveColumns.clear();
				for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
				{
					if (psolutionMethod_->columns_[i]->isInCurrentSolution())
					{
						positiveColumns.push_back(psolutionMethod_->columns_[i]);
					}
				}

				if (has_compatible_column)
				{
					// On recalcule les degr�s d'incompatibilit�s
					IncompatibilityDegree id(psolutionMethod_, positiveColumns, psolutionMethod_->tasks_);
					std::vector<IB_Column *> columns_to_recompute;
					for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
					{
						IB_Column *column = psolutionMethod_->columns_[i];
						if (!column->isInCurrentSolution())
						{
							bool recompute = column->getPhase() == -1;
							if (recompute)
							{
								columns_to_recompute.push_back(column);
							}
						}
					}

					int n_threads = 8;

					if (columns_to_recompute.size() <= 100)
					{
						n_threads = 1;
					}

					std::vector<std::vector<IB_Column *>> thread_tasks;
					std::vector<IB_Column *> current_vector;
					for (int i = 0; i < columns_to_recompute.size(); i++)
					{
						if (current_vector.size() >= ((float)columns_to_recompute.size()) / n_threads)
						{
							thread_tasks.push_back(current_vector);
							current_vector = {};
						}

						current_vector.push_back(columns_to_recompute[i]);
					}

					if (current_vector.size())
					{
						thread_tasks.push_back(current_vector);
					}

					std::cout << thread_tasks.size() << " est la taille." << std::endl;

					std::vector<std::thread> threads_;
					for (auto thread_task : thread_tasks)
					{
						threads_.push_back(std::thread(calcIncompatibilityDegrees, id, thread_task));
					}

					for (int i = 0; i < threads_.size(); i++)
					{
						threads_[i].join();
					}
				}

				continue_ = getBCompatibleColumnId(&colsIn, &colsOut);
			}
		}

		colsIn.clear();
		colsOut.clear();
		// R�solution du probl�me compl�mentaire
		std::vector<int> phaseSeq = getPhaseSequence();

		bool integral = false;
		std::vector<double> solution;

		double pastCost = currentCost_;
		std::vector<double> duals;
		for (int i = 0; i < phaseSeq.size(); i++)
		{
			if (skipPhase && (phaseSeq[i] != -1 && phaseSeq[i] <= previousPhaseMax))
			{
				continue;
			}
			phaseMax = phaseSeq[i];
			std::cout << "Resolution du probleme complementaire en phase : " << phaseSeq[i] << std::endl;
			solution.clear();
			
			IB_ComplementaryProblem cp(psolutionMethod_, phaseSeq[0]);
			cp.setPhase(phaseSeq[i]);
			cp.constructProblem();
			duals.clear();
			double objective = cp.solve(&solution, NULL, &duals);
			double pas = 100000;
			for(int kp = 0; kp < solution.size();kp++) {
				if(fabs(solution[kp]) > 1e-5 && 1./solution[kp] < pas) {
					pas = 1./solution[kp];

				}
			}
			
			last_objective = objective;
			double gapValue = fabs(objective * psolutionMethod_->sum_bi) / currentCost_;
			double gapValue2 = fabs(objective * pas) / (currentCost_ - bound);
			std::cout << gapValue << std::endl;
			if ((objective >= 0 || gapValue <= 0.005) && phaseSeq[i] == -1)
			{
				std::cout << "Echec du probleme complementaire en phase -1" << std::endl;
				solved = true;
				cp.destroy();
				return;
			}

			if (phaseSeq[i] != -1 && (objective >= 0 || gapValue2 <= 0.006))
			{
				std::cout << "Echec du probleme complementaire" << std::endl;

				cp.destroy();
				continue;
			}

			if (gapValue2 >= 0.006 && isIntegral(solution))
			{
				colsIn.clear();
				colsOut.clear();
				integral = true;
				while (objective <= 0 && gapValue2 >= 0.006 && isIntegral(solution))
				{
					std::cout << "Direction entiere trouvee " << std::endl;

					for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
					{
						IB_Column *column = psolutionMethod_->columns_[i];
						if (solution[i] > 1e-5)
						{
							if (column->isInCurrentSolution())
							{
								colsOut.push_back(i);
							}
							else
							{
								colsIn.push_back(i);
							}
						}
					}

					solution.clear();
					objective = cp.solve(&solution);
					gapValue = fabs(objective * psolutionMethod_->sum_bi) / currentCost_;
					pas = 100000;
					for(int kp = 0; kp < solution.size();kp++) {
						if(fabs(solution[kp]) > 1e-5 && 1./solution[kp] < pas) {
							pas = 1./solution[kp];

						}
					}

					gapValue2 = fabs(objective * pas) / (currentCost_ - bound);
				}

				if (integral)
				{
					pivotColumnsInSolution(colsIn, colsOut);
				}
				cp.destroy();
				break;
			}

			if (gapValue2 >= 0.006)
			{
				if (!integral)
				{
					cp.destroy();
				}
				artificial_column_id += 1;
				IB_Column *artificial_column = new IB_Column("artificial_column_" + std::to_string(artificial_column_id),
															 std::vector<std::pair<std::string, int>>());
				std::vector<int> acolsOut, acolsIn;
				for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
				{
					IB_Column *column = psolutionMethod_->columns_[i];
					if (solution[i] > 1e-5)
					{
						if (column->isInCurrentSolution())
						{
							acolsOut.push_back(i);
						}
						else
						{
							acolsIn.push_back(i);
						}
					}
				}

				if (addColumns_ && canBeInIntegerDirection(artificial_column, solution))
				{
					n_added_columns += 1;
					columnAdded = true;
					int colId = psolutionMethod_->addArtificialColumn(artificial_column);
					colsIn.clear();
					colsOut.clear();
					auto add_columns_start = std::chrono::high_resolution_clock::now();
					std::pair<bool, int> result = cpWithArtificialColumn(colId, std::set<int>(acolsIn.begin(), acolsIn.end()),
																		 std::set<int>(acolsOut.begin(), acolsOut.end()), colsIn, colsOut, phaseSeq[i], 3, objective / (acolsIn.size() + acolsOut.size()));
					int add_columns_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - add_columns_start).count();
					n_added_columns_it = result.second;
					if (result.first)
					{
						double reduced_cost_rc = 0;
						for (auto col : colsIn)
						{
							reduced_cost_rc += psolutionMethod_->columns_[col]->getCost();
						}

						for (auto col : colsOut)
						{
							reduced_cost_rc -= psolutionMethod_->columns_[col]->getCost();
						}

						if (compete)
						{
							std::vector<int> newColsIn, newColsOut;
							auto zoom_start = std::chrono::high_resolution_clock::now();
							bool zoomSuccess = zoom(phaseMax, phaseSeq, solution, &newColsIn, &newColsOut);
							if (zoomSuccess)
							{
								double reduced_cost = 0;
								for (auto col : newColsIn)
								{
									reduced_cost += psolutionMethod_->columns_[col]->getCost();
								}

								for (auto col : newColsOut)
								{
									reduced_cost -= psolutionMethod_->columns_[col]->getCost();
								}

								int zoom_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - zoom_start).count();

								addCompeteRow(reduced_cost_rc, add_columns_time, reduced_cost, zoom_time);
							}
						}
						std::cout << "bon c'est bon" << std::endl;
						IB_CompatibilityChecker compatibilityChecker(psolutionMethod_);
						std::vector<int> nonNullColumnsIndices;
						for (auto colInd : colsIn)
						{
							nonNullColumnsIndices.push_back(colInd);
						}

						for (auto colInd : colsOut)
						{
							nonNullColumnsIndices.push_back(colInd);
						}
						pivot_distance = compatibilityChecker.getDistance(nonNullColumnsIndices);
						// psolutionMethod_->removeColumn(colId);
						delete artificial_column;
						artificial_column = NULL;
						integral = true;
						n_success_add_columns += 1;
						pivotColumnsInSolution(colsIn, colsOut);
						break;
					} else {
						if(compete) {
							addCompeteRow(0, add_columns_time, 0, 0);
						}
					}

					// psolutionMethod_->removeColumn(colId);
					delete artificial_column;
					artificial_column = NULL;
					break;
				}
				else
				{
					delete artificial_column;
					artificial_column = NULL;
					break;
				}
			}
			else if (integral)
			{
				break;
			}
		}
		skipPhase = false;

		if (!integral)
		{
			bool do_zoom = true;
			if (disEnabled)
			{
				std::vector<int> acolsIn, acolsOut;
				std::map<int, int> colsCorrespondance;
				std::map<int, int> colsCorrespondanceInv;
				std::vector<IB_Column *> involvedColumns;
				for (int i = 0; i < psolutionMethod_->columns_.size(); i++)
				{
					IB_Column *column = psolutionMethod_->columns_[i];
					if (column->isInCurrentSolution() || (phaseMax == -1 || column->getPhase() <= phaseMax))
					{
						involvedColumns.push_back(column);
						colsCorrespondance[i] = involvedColumns.size() - 1;
						colsCorrespondanceInv[involvedColumns.size() - 1] = i;
					}
				}

				for (int i = 0; i < solution.size(); i++)
				{
					if (solution[i] > 1e-5)
					{
						if (psolutionMethod_->columns_[i]->isInCurrentSolution())
						{
							acolsOut.push_back(colsCorrespondance[i]);
						}
						else
						{
							acolsIn.push_back(colsCorrespondance[i]);
						}
					}
				}
				colsIn.clear();
				colsOut.clear();
				ISUD_Base problem(psolutionMethod_->tasks_, psolutionMethod_->rhs_, involvedColumns);
				std::pair<bool, int> cpDisResult = std::pair<bool, int>(false, 1);
				int add_columns_time;
				double reduced_cost_rc = 0;
				int disp_size = 0;
				double remaining = currentCost_ - bound;
				if (phaseMax != -1)
				{
					dis_problem_size = "";
					auto disaggregation_start = std::chrono::high_resolution_clock::now();
					cpDisResult = cpWithDisaggregation(duals, phaseMax, &colsIn, &colsOut, acolsOut, acolsIn, &problem, colsCorrespondanceInv, last_objective, psolutionMethod_->columns_.size(), &disp_size, bound);
					
					add_columns_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - disaggregation_start).count();

				}
				if (!cpDisResult.first && cpDisResult.second == 2)
				{
					columnAdded = true;
					n_added_columns += 1;
				}
				if (cpDisResult.first)
				{
					n_added_columns += 1;
					if (compete)
					{
						for (auto col : colsIn)
						{
							reduced_cost_rc += psolutionMethod_->columns_[col]->getCost();
						}

						for (auto col : colsOut)
						{
							reduced_cost_rc -= psolutionMethod_->columns_[col]->getCost();
						}

						std::vector<int> newColsIn, newColsOut;
						auto zoom_start = std::chrono::high_resolution_clock::now();
						bool zoomSuccess = zoom(phaseMax, phaseSeq, solution, &newColsIn, &newColsOut);
						if (zoomSuccess)
						{
							double reduced_cost = 0;
							for (auto col : newColsIn)
							{
								reduced_cost += psolutionMethod_->columns_[col]->getCost();
							}

							for (auto col : newColsOut)
							{
								reduced_cost -= psolutionMethod_->columns_[col]->getCost();
							}

							int zoom_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - zoom_start).count();

							addCompeteRow(reduced_cost_rc, add_columns_time, reduced_cost, zoom_time, disp_size, last_objective, remaining);
						}
					}

					IB_CompatibilityChecker compatibilityChecker(psolutionMethod_);
					std::vector<int> nonNullColumnsIndices;
					for (auto colInd : colsIn)
					{
						nonNullColumnsIndices.push_back(colInd);
					}

					for (auto colInd : colsOut)
					{
						nonNullColumnsIndices.push_back(colInd);
					}
					pivot_distance = 0;

					pivotColumnsInSolution(colsIn, colsOut);
					columnAdded = true;
					do_zoom = false;
					n_success_add_columns += 1;
				} else if(compete && cpDisResult.second == 2) {
					addCompeteRow(0, add_columns_time, 0, 0, disp_size, last_objective, remaining);
				}
			}
			// solution.clear();
			/*IB_ComplementaryProblem cp2(psolutionMethod_, phaseSeq[0]);
			if (phaseMax != -1 and phaseMax < 5) {
				cp2.setPhase(5);
			}
			double objective = cp2.solve(&solution);*/
			if (do_zoom)
			{
				colsIn.clear();
				colsOut.clear();
				// On r�cup�re une meilleure solution gr�ce � zoom
				bool zoomSuccess = zoom(phaseMax, phaseSeq, solution, &colsIn, &colsOut);
				has_zoom = true;
				if (zoomSuccess)
				{
					pivotColumnsInSolution(colsIn, colsOut);
				}
				else if (phaseMax == -1)
				{
					solved = true;
					return;
				}
				else
				{
					skipPhase = true;
					previousPhaseMax = phaseMax;
				}
			}
		}

		double newCost = currentCost_;
		n_iterations += 1;
		auto current_time = std::chrono::high_resolution_clock::now();
		addLine(newCost, newCost - pastCost, n_iterations, n_added_columns, n_success_add_columns, has_zoom, phaseMax, columnAdded, n_added_columns_it, pivot_distance,
				std::chrono::duration_cast<std::chrono::milliseconds>(current_time - iteration_start).count(),
				std::chrono::duration_cast<std::chrono::milliseconds>(current_time - global_start).count());

		if ((newCost - bound) / newCost <= 0.005)
		{
			return;
		}
	}
}

// Complementary problem with task disaggregation
std::pair<bool, int> ISUD::cpWithDisaggregation(std::vector<double>& duals, int phase, std::vector<int> *colsIn, std::vector<int> *colsOut, std::vector<int> &acolsOut, std::vector<int> &acolsIn, ISUD_Base *problem, std::map<int, int> originalProblemColumns, double past_objective,
												int n_cols, int* size_dis_problem, double bound)
{
	DisComputer dis_computer(problem, originalProblemColumns, phase);
	std::map<int, std::map<int, int>> posColsAffectation;
	double objective;
	std::map<int, std::vector<int>> disaggregation = dis_computer.computeDisaggregation(acolsIn, acolsOut, &posColsAffectation, &objective);

	std::cout << "Desagregation" << std::endl;
	int n_dis = 0;
	bool has_dis = false;
	for (auto pair : disaggregation)
	{
		if (pair.second.size() > 1)
		{
			n_dis += 1;
			std::cout << "Tache " << pair.first << std::endl;
			for (auto rhs : pair.second)
			{
				std::cout << "     RHS " << rhs << std::endl;
			}

			has_dis = true;
		}
	}
	std::cout << "Objectif : " << objective << " " << past_objective << std::endl;
	if (!has_dis)
	{
		return std::pair<bool, int>(false, 1);
	}

	std::vector<int> involvedColumns;
	for (int i = 0; i < problem->columns_.size(); i++)
	{
		involvedColumns.push_back(i);
	}

	DisProblem problem_dis(problem, involvedColumns, disaggregation);
	int n_zero_cols = problem_dis.constructColumns(originalProblemColumns, phase, &posColsAffectation);
	*size_dis_problem = n_zero_cols;
	std::vector<IB_Column *> new_columns = problem_dis.getAllColumns();
	dis_problem_size = std::to_string(problem->tasks_.size()) + ", " + std::to_string(problem->columns_.size()) + " | " +
					   std::to_string(problem_dis.newTasks.size()) + ", " + std::to_string(new_columns.size());
	if(*size_dis_problem >= 15000) {
		return std::pair<bool, int>(false, 1);
	}
	// On calcule les degrés d'incompatibilités et on appelle le problème complémentaire
	/*
	std::vector<IB_Column*> positiveColumns;
	for (auto column : new_columns) {
		if (column->isInCurrentSolution()) {
			positiveColumns.push_back(column);
		}
	}

	ISUD_Base* np = problem_dis.getISUDBase();
	std::cout << new_columns.size() << " colonnes dans le probleme complementaire desagrege" << std::endl;
	std::cout << "Calcul des degrés d'incompatibilités" << std::endl;
	IncompatibilityDegree id(np, positiveColumns, np->tasks_);
	for (auto column : new_columns) {
		column->setPhase(id.getIncompatibilityDegree(column));
	}*/

	ISUD_Base *np = problem_dis.getISUDBase();

	std::cout << "Obtention des colonnes obligatoires" << std::endl;

	std::vector<int> compulsoryCols;
	ComplementaryProblemDis cp(np, -1);
	cp.constructProblem();
	std::vector<double> solution;
	double final_objective = cp.solve(&solution, problem_dis.originalProblemColumns_, problem_dis.tasksMapping, duals, true);
	std::cout << "Objectif final : " << final_objective << std::endl;
	std::vector<double> solutionAggregated = problem_dis.aggregate(solution, n_cols);
	double initial_objective = final_objective;

	if (isIntegral(solutionAggregated) && final_objective < -1e-4)
	{

		while (isIntegral(solutionAggregated) && final_objective < -1e-4)
		{
			std::cout << "Direction entiere trouvee " << std::endl;

			for (int i = 0; i < solutionAggregated.size(); i++)
			{
				if (solutionAggregated[i] > 1e-5)
				{
					if (psolutionMethod_->columns_[i]->isInCurrentSolution())
					{
						colsOut->push_back(i);
					}
					else
					{
						colsIn->push_back(i);
					}
				}
			}

			solution.clear();
			final_objective = cp.solve(&solution, problem_dis.originalProblemColumns_, problem_dis.tasksMapping, duals, false);
			
			solutionAggregated = problem_dis.aggregate(solution, n_cols);
		}

		std::cout << "DESAGREGATION REUSSIE" << std::endl;
		problem_dis.deleteProblem();
		return std::pair<bool, int>(true, 0);
	}

	std::vector<int> ncolsOut, ncolsIn;
	for (int i = 0; i < solution.size(); i++)
	{
		if (solution[i] > 1e-5)
		{
			if (new_columns[i]->isInCurrentSolution())
			{
				ncolsOut.push_back(i);
			}
			else
			{
				ncolsIn.push_back(i);
			}
		}
	}

	problem_dis.deleteProblem();

	return std::pair<bool, int>(false, 2);
	;

	return cpWithDisaggregation(duals, phase, colsIn, colsOut, ncolsOut, ncolsIn, np, problem_dis.originalProblemColumns_, final_objective, n_cols, size_dis_problem, bound);
}


// Return phase sequence for multiphase strategy
std::vector<int> ISUD::getPhaseSequence()
{
	std::vector<int> incompatibilityDegrees;
	for (auto column : psolutionMethod_->columns_)
	{
		if (!column->isInCurrentSolution())
		{
			incompatibilityDegrees.push_back(column->getPhase());
		}
	}
	std::sort(incompatibilityDegrees.begin(), incompatibilityDegrees.end());
	std::vector<int> phaseSeq = {incompatibilityDegrees[int(floor(incompatibilityDegrees.size() * 0.05))],
								 incompatibilityDegrees[int(floor(incompatibilityDegrees.size() * 0.1))],
								 incompatibilityDegrees[int(floor(incompatibilityDegrees.size() * 0.3))],
								 incompatibilityDegrees[int(floor(incompatibilityDegrees.size() * 0.4))],
								 incompatibilityDegrees[int(floor(incompatibilityDegrees.size() * 0.5))],
								 incompatibilityDegrees[int(floor(incompatibilityDegrees.size() * 0.6))],
								 -1};

	if (phaseSeq[0] == 0)
	{
		phaseSeq[0] = 1;
	}

	std::vector<int> finalPhaseSeq;
	for (int i = 0; i < phaseSeq.size(); i++)
	{
		if (std::find(finalPhaseSeq.begin(), finalPhaseSeq.end(), phaseSeq[i]) == finalPhaseSeq.end())
		{
			finalPhaseSeq.push_back(phaseSeq[i]);
		}
	}

	// finalPhaseSeq = std::vector<int>({ 1, 2, 3, 4, 5, 6, -1 });
	// return {-1};
	return finalPhaseSeq;
	// return std::vector<int>({ 1, 2, 3, 4, 5, 6, -1 });
}
