#include "DisProblem.h"
#include <string>
#include "IncompatibilityDegree.h"
#include <thread>
#include <future>


void computeColumns(int phase, std::vector<IB_Column*> columns, ISUD_Base* intermediateProblem, std::vector<IB_Column*> positiveColumns, std::vector<std::string> newTasks,
	std::map<std::pair<std::string, std::string>, std::pair<int, int>> tasksCorrespondance2, std::promise<std::vector< std::vector<std::pair<std::map<int, int>, int>>>>&& promiseReturnValue) {
	IncompatibilityDegree id(intermediateProblem, positiveColumns, newTasks, &tasksCorrespondance2);

	std::vector< std::vector<std::pair<std::map<int, int>, int>>> returnValue;

	for (auto contributions : columns) {
		returnValue.push_back(id.getAllPaths(contributions, phase));
	}

	promiseReturnValue.set_value(returnValue);
}

int DisProblem::constructColumns(std::map<int, int> originalProblemColumns, int phase, std::map<int, std::map<int, int>>* posColumnsAffectation) {
	std::vector<IB_Column*>& initialColumns = psolutionMethod_->columns_;
	allColumns.empty();

	std::map<std::pair<int, int>, int> remainingRhs;
	std::set<int> disaggregatedTasks;
	int n_tasks = 0;
	for (auto pair : rhs_mapping) {
		int initial_task = pair.first;
		std::vector<int> mappedRhs = pair.second;
		int j = 0;
		if(mappedRhs.size() > 1) {
			disaggregatedTasks.insert(initial_task);
		}
		for (auto rhs : mappedRhs) {
			remainingRhs[std::pair<int, int>(initial_task, j)] = rhs;
			std::string taskName = "task_" + std::to_string(n_tasks);
			newRhs.push_back(rhs);
			n_tasks += 1;
			tasksMapping[taskName] = std::pair<int, int>(initial_task, j);
			tasksMappingInv[std::pair<int, int>(initial_task, j)] = taskName;
			newTasks.push_back(taskName);
			j += 1;
		}
	}

	


	int col_nb = 0;
	std::vector<int> zeroColumns;
	std::vector<int> otherColumns;
	std::vector<IB_Column*> positiveColumns;
	std::map<std::pair<std::string, std::string>, std::pair<int, int>> tasksCorrespondance2;
	for (int i = 0; i < involvedColumns_.size(); i++) {
		int indice = involvedColumns_[i];
		IB_Column* column = initialColumns[indice];
		if (column->isInCurrentSolution()) {
			positiveColumns.push_back(column);
			std::vector<std::map<int, int>> tasksInstances = getTasksInstances(indice, column, &remainingRhs, NULL);
			for (auto taskInstance : tasksInstances) {
				std::vector<std::pair<std::string, int>> contribs;

				for (auto contrib : column->getContribs()) {
					contribs.push_back(std::pair<std::string, int>(tasksMappingInv[std::pair<int, int>(tasksCorrespondance[contrib.first], taskInstance[tasksCorrespondance[contrib.first]])], contrib.second));
					tasksCorrespondance2[std::pair<std::string, std::string>(contrib.first, column->getName())] = std::pair<int, int>(tasksCorrespondance[contrib.first], taskInstance[tasksCorrespondance[contrib.first]]);
				}

				IB_Column* new_column = new IB_Column("column_" + std::to_string(col_nb), contribs);
				new_column->setCost(column->getCost());
				new_column->reorderContribs(newTasks, true);
				new_column->setInCurrentSolution();
				new_column->setPhase(0);

				allColumns.push_back(new_column);
				columnsCorrespondance[allColumns.size() - 1] = taskInstance;
				initialColumnAssociated[allColumns.size() - 1] = indice;
				originalProblemColumns_[allColumns.size() - 1] = originalProblemColumns[indice];
				if (initialColumnAssociatedInv.find(indice) == initialColumnAssociatedInv.end()) {
					initialColumnAssociatedInv[indice] = std::vector<int>();
				}

				initialColumnAssociatedInv[indice].push_back(allColumns.size() - 1);
				col_nb += 1;
			}
		}
		else {
			bool take = false;
			for(auto contrib: column->getContribs()) {
				if(contrib.second != 0) {
					if(disaggregatedTasks.count(tasksCorrespondance[contrib.first])) {
						take = true;
					}
				}
			}
			if(take) {
				zeroColumns.push_back(indice);
			} else {
				otherColumns.push_back(indice);
			}
		}
	}

	if(zeroColumns.size() >= 15000) {
		return zeroColumns.size();
	}

	int n_threads = 8;
	std::vector<std::vector<IB_Column*>> threads_tasks;
	std::vector<std::vector<int>> threads_indices;
	std::vector<IB_Column*> current_vector;
	std::vector<int> current_vector_indices;
	std::vector<std::promise<std::vector< std::vector<std::pair<std::map<int, int>, int>>>>> return_values;
	for (int i = 0; i < zeroColumns.size(); i++) {
		if (current_vector.size() >= ((float)zeroColumns.size()) / n_threads && current_vector.size() > 0) {
			threads_tasks.push_back(current_vector);
			return_values.push_back(std::promise<std::vector< std::vector<std::pair<std::map<int, int>, int>>>>());
			threads_indices.push_back(current_vector_indices);
			current_vector = {};
			current_vector_indices = {};
		}

		int indice = zeroColumns[i];
		IB_Column* column = initialColumns[indice];
		
		current_vector.push_back(column);
		current_vector_indices.push_back(indice);
	}

	if (current_vector.size()) {
		threads_tasks.push_back(current_vector);
		return_values.push_back(std::promise<std::vector< std::vector<std::pair<std::map<int, int>, int>>>>());
		threads_indices.push_back(current_vector_indices);
	}


	ISUD_Base intermediateProblem(newTasks, newRhs, allColumns);
	
	std::vector<std::thread> threads_;
	std::vector< std::future<std::vector<std::vector<std::pair<std::map<int, int>, int>>>>> rv_futures;
	for (int i = 0; i < threads_tasks.size(); i++) {
		rv_futures.push_back(return_values[i].get_future());

		threads_.push_back(std::thread(&computeColumns, phase, threads_tasks[i], psolutionMethod_, positiveColumns, psolutionMethod_->tasks_, tasksCorrespondance2, std::move(return_values[i])));
	}

	for (int i = 0; i < threads_.size(); i++) {
		threads_[i].join();
	}

	std::vector<int> zeroColumnsIndices;
	std::vector< std::vector<std::pair<std::map<int, int>, int>> > results;
	for (int i = 0; i < return_values.size(); i++) {
		std::vector< std::vector<std::pair<std::map<int, int>, int>> > newResults = rv_futures[i].get();
		results.insert(results.end(), newResults.begin(), newResults.end());
		zeroColumnsIndices.insert(zeroColumnsIndices.end(), threads_indices[i].begin(), threads_indices[i].end());
	}

	//IncompatibilityDegree id(&intermediateProblem, allColumns, newTasks);
	for (int i = 0; i < zeroColumnsIndices.size(); i++) {
		int indice = zeroColumnsIndices[i];
		IB_Column* column = initialColumns[indice];
		if (!column->isInCurrentSolution()) {
			/*std::set<int> contributions__;
			for (auto contrib : column->getContribs()) {
				if (contrib.second != 0) {
					contributions__.insert(tasksCorrespondance[contrib.first]);
				}
			}*/
			 
			std::vector<std::pair<std::map<int, int>, int>> tasksInstances = results[i];
			if(tasksInstances.size() == 0) {
				std::cout << "La phase etait : " << column->getPhase() << std::endl;
			}
			for (auto taskInstancePair : tasksInstances) {
				auto taskInstance = taskInstancePair.first;

				std::vector<std::pair<std::string, int>> contribs;

				for (auto contrib : column->getContribs()) {
					contribs.push_back(std::pair<std::string, int>(tasksMappingInv[std::pair<int, int>(tasksCorrespondance[contrib.first], taskInstance[tasksCorrespondance[contrib.first]])], contrib.second));
				}

				IB_Column* new_column = new IB_Column("column_" + std::to_string(col_nb), contribs);
				new_column->setCost(column->getCost());
				new_column->reorderContribs(newTasks, true);
				new_column->setOutCurrentSolution();
				if (taskInstancePair.second == -1) {
					new_column->setPhase(column->getPhase());
				} else{
                       new_column->setPhase(taskInstancePair.second);
				}

				allColumns.push_back(new_column);
				columnsCorrespondance[allColumns.size() - 1] = taskInstance;
				initialColumnAssociated[allColumns.size() - 1] = indice;
				originalProblemColumns_[allColumns.size() - 1] = originalProblemColumns[indice];
				if (initialColumnAssociatedInv.find(indice) == initialColumnAssociatedInv.end()) {
					initialColumnAssociatedInv[indice] = std::vector<int>();
				}

				initialColumnAssociatedInv[indice].push_back(allColumns.size() - 1);
				col_nb += 1;
			}
		}
	}

	//IncompatibilityDegree id(&intermediateProblem, allColumns, newTasks);
	for (int i = 0; i < otherColumns.size(); i++) {
		int indice = otherColumns[i];
		IB_Column* column = initialColumns[indice];
		if (!column->isInCurrentSolution()) {
			/*std::set<int> contributions__;
			for (auto contrib : column->getContribs()) {
				if (contrib.second != 0) {
					contributions__.insert(tasksCorrespondance[contrib.first]);
				}
			}*/
			 
			std::map<int, int> taskInstance;
			for(auto contrib: column->getContribs()) {
				if(contrib.second != 0) taskInstance[tasksCorrespondance[contrib.first]] = 0;
			}
			std::vector<std::pair<std::map<int, int>, int>> tasksInstances = { std::pair<std::map<int, int>, int>(taskInstance, column->getPhase())};
			
			for (auto taskInstancePair : tasksInstances) {
				auto taskInstance = taskInstancePair.first;

				std::vector<std::pair<std::string, int>> contribs;

				for (auto contrib : column->getContribs()) {
					contribs.push_back(std::pair<std::string, int>(tasksMappingInv[std::pair<int, int>(tasksCorrespondance[contrib.first], taskInstance[tasksCorrespondance[contrib.first]])], contrib.second));
				}

				IB_Column* new_column = new IB_Column("column_" + std::to_string(col_nb), contribs);
				new_column->setCost(column->getCost());
				new_column->reorderContribs(newTasks, true);
				new_column->setOutCurrentSolution();
				if (taskInstancePair.second == -1) {
					new_column->setPhase(column->getPhase());
				} else{
                       new_column->setPhase(taskInstancePair.second);
				}

				allColumns.push_back(new_column);
				columnsCorrespondance[allColumns.size() - 1] = taskInstance;
				initialColumnAssociated[allColumns.size() - 1] = indice;
				originalProblemColumns_[allColumns.size() - 1] = originalProblemColumns[indice];
				if (initialColumnAssociatedInv.find(indice) == initialColumnAssociatedInv.end()) {
					initialColumnAssociatedInv[indice] = std::vector<int>();
				}

				initialColumnAssociatedInv[indice].push_back(allColumns.size() - 1);
				col_nb += 1;
			}
		}
	}


	/*
	int col_nb = 0;
	for (int i = 0; i < involvedColumns_.size(); i++) {
		int indice = involvedColumns_[i];
		IB_Column* column = initialColumns[indice];
		std::vector<std::map<int, int>> tasksInstances = getTasksInstances(indice, column, &remainingRhs, posColumnsAffectation);
		for (auto taskInstance : tasksInstances) {
			std::vector<std::pair<std::string, int>> contribs;

			for (auto contrib : column->getContribs()) {
				contribs.push_back(std::pair<std::string, int>(tasksMappingInv[std::pair<int, int>(tasksCorrespondance[contrib.first], taskInstance[tasksCorrespondance[contrib.first]])], contrib.second));
			}

			IB_Column* new_column = new IB_Column("column_" + std::to_string(col_nb), contribs);
			new_column->setCost(column->getCost());
			new_column->reorderContribs(newTasks);
			if (column->isInCurrentSolution()) {
				new_column->setInCurrentSolution();
			}
			else {
				new_column->setOutCurrentSolution();
			}

			allColumns.push_back(new_column);
			columnsCorrespondance[allColumns.size() - 1] = taskInstance;
			initialColumnAssociated[allColumns.size() - 1] = indice;
			originalProblemColumns_[allColumns.size() - 1] = originalProblemColumns[indice];
			if (initialColumnAssociatedInv.find(indice) == initialColumnAssociatedInv.end()) {
				initialColumnAssociatedInv[indice] = std::vector<int>();
			}

			initialColumnAssociatedInv[indice].push_back(allColumns.size() - 1);
			col_nb += 1;
		}
	}*/

	current_problem = new ISUD_Base(newTasks, newRhs, allColumns);
	return zeroColumns.size();
}

ISUD_Base* DisProblem::getISUDBase() {
	return current_problem;
}

std::vector<int> DisProblem::getAssociatedColumns(int task, int n, int column_no) {
	std::vector<int> returnValue;
	int indice = initialColumnAssociated[column_no];
	std::vector<int> associatedColumns = initialColumnAssociatedInv[indice];
	for (int colId : associatedColumns) {
		for (auto taskInstance : columnsCorrespondance[colId]) {
			if (taskInstance.first == task && taskInstance.second == n) {
				returnValue.push_back(colId);
			}
		}
	}

	return returnValue;
}

std::pair<int, int> DisProblem::getPosition(std::string taskName) {
	return tasksMapping[taskName];
}

std::vector<std::map<int, int>> DisProblem::getTasksInstances(int colId, IB_Column* column, std::map<std::pair<int, int>, int>* remainingRhs, std::map<int, std::map<int, int>>* posColumnsAffectation) {
	if (column->isInCurrentSolution()) {
		// Si la colonne est dans la solution courante, on garde une seule colonne
		std::map<int, int> correspondance;
		
		for (auto contrib : column->getContribs()) {
				int taskId = tasksCorrespondance[contrib.first];
				int n_max = rhs_mapping[taskId].size();

				bool found = false;
				for (int j = 0; j < n_max; j++) {
					if ((*remainingRhs)[std::pair<int, int>(taskId, j)] > 0) {
						int new_val = (*remainingRhs)[std::pair<int, int>(taskId, j)] - 1;
						remainingRhs->erase(std::pair<int, int>(taskId, j));
						remainingRhs->emplace(std::pair<int, int>(taskId, j), new_val);
						correspondance[taskId] = j;
						break;
					}
				}

			}

		return std::vector<std::map<int, int>>({ correspondance });
	}
	else {
		std::vector<std::pair<std::string, int>> contribs;
		for (auto contrib : column->getContribs()) {
			contribs.push_back(std::pair<std::string, int>(contrib.first, contrib.second));
		}

		




		return getTasksInstancesBis(contribs);
	}
}

std::vector<std::map<int, int>> DisProblem::getTasksInstancesBis(std::vector<std::pair<std::string, int>> contribs) {
	std::vector<std::map<int, int>> returnValue;
	std::pair<std::string, int> contrib = contribs[0];
	contribs.erase(contribs.begin());
	for (int i = 0; i < rhs_mapping[tasksCorrespondance[contrib.first]].size(); i++) {
		if (contribs.size() > 0) {
			for (auto corr : getTasksInstancesBis(contribs)) {
				std::map<int, int> correspondance;
				correspondance[tasksCorrespondance[contrib.first]] = i;
				for (auto pair : corr) {
					correspondance[pair.first] = pair.second;
				}

				returnValue.push_back(correspondance);
			}
		}
		else {
			std::map<int, int> correspondance;
			correspondance[tasksCorrespondance[contrib.first]] = i;
			returnValue.push_back(correspondance);
		}
		
	}

	return returnValue;
}

std::vector<double> DisProblem::aggregate(std::vector<double>& solution, int n_cols) {
	std::vector<double> return_value;
	for (int i = 0; i < n_cols; i++) {
		return_value.push_back(0);
	}

	for (int i = 0; i < solution.size(); i++) {
		return_value[originalProblemColumns_[i]] += solution[i];
	}

	return return_value;
}

void DisProblem::deleteProblem() {
	current_problem->destroyColumns();
	delete current_problem;
	current_problem = NULL;
	allColumns.empty();
}