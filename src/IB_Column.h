#ifndef IB_Column_H
#define IB_Column_H
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <set>

// Represent a column of ISUD
class IB_Column {
    private:
        // If the column is compatible, indices of positive columns making this column compatible
        std::vector<std::string> compatibleBy_;
        // Contributions of the columns (task to number)
        std::vector<std::pair<std::string, int>> contributions;
        // Contributions map
        std::map<std::string, int> contribsIndex_;
        // Tasks covered by column (set)
        std::set<std::string> tasksInColumn;
        // Column state (for example column compatible)
        int column_state_;
        // Is the column positive ?
        bool is_in_current_solution_ = false;
        // Cost of the column
        double cost_;
        // Inc degree of the column
        int phase_;
        // Name of the column
        std::string name_;
        // If column is positive
        bool is_in_p = false;
        // Recompute the state of the column
        bool recompute_compatibility = true;
        // Order of tasks covered by the column
        std::vector<std::string> contribsOrder;
        // Is column positive and in the set of linearly independant columns ?
        bool is_in_p_independent = false;

    public:
        // STATES of column
        static const int STATE_COMPATIBLE = 1;
        static const int STATE_INCOMPATIBLE = 2;
        static const int STATE_COMPATIBLE_LINEARLY = 3;

        // Constructor of column
        IB_Column (std::string name, std::vector<std::pair<std::string, int>> contribs) : name_(name), contributions(contribs), cost_(0), compatibleBy_() {
            compatibleBy_.clear();
            column_state_ = STATE_INCOMPATIBLE;
            is_in_current_solution_ = false;

            for (auto pair : contribs) {
                contribsIndex_[pair.first] = pair.second;
                tasksInColumn.insert(pair.first);
            }
        };

        // Are columns covering exactly the tasks coveredTasks
        bool isSame(std::set<std::string> coveredTasks) {
            return tasksInColumn == coveredTasks;
        }

        // Set compatible by columns in cb
        void setCompatibleby(std::set<std::string>  cb) {
            compatibleBy_.clear();
            
            for (auto colIndice : cb) {
                compatibleBy_.push_back(colIndice);
            }
        }

        // Reorder contributions of columns with task_order
        void reorderContribs(std::vector<std::string> tasks_order, bool same = false) {
            contribsOrder.clear();

            if (!same) {
                for (int i = 0; i < tasks_order.size(); i++) {
                    if (contribsIndex_.find(tasks_order[i]) != contribsIndex_.end()) {
                        contribsOrder.push_back(tasks_order[i]);
                    }
                }
            }
            else {
                for (auto contrib : contributions) {
                    if (contrib.second != 0) {
                        contribsOrder.push_back(contrib.first);
                    }
                }
            }
            

        }

        // Get columns compatible by
        std::vector<std::string> getCompatibleBy() {
            return compatibleBy_;
        }

        // Column name
        std::string getName() {
            return name_;
        }
        
        // Set column name
        void setName(std::string name) {
            name_ = name;
        }

        // Change state
        void setCompatible() {
            column_state_ = STATE_COMPATIBLE;
        }

        void setIncompatible() {
            column_state_ = STATE_INCOMPATIBLE;
        }

        void setCompatibleLinearly() {
            column_state_ = STATE_COMPATIBLE_LINEARLY;
        }

        // Put in or out current solution
        void setInCurrentSolution() {
            is_in_current_solution_ = true;
        }

        void setOutCurrentSolution() {
            is_in_current_solution_ = false;
        }

        // Put column positive
        void changeInP(bool newState) {
            is_in_p = newState;
        }

        // Is coluumn positive ?
        bool isInP() {
            return is_in_p;
        }

        // Put column in positive independent columns set
        void changeInPIndependent(bool new_state) {
            is_in_p_independent = new_state;
        }

        bool isInPIndependent() {
            return is_in_p_independent;
        }

        // Do we need to recompute the state ?
        bool shouldRecomputeCompatibility() {
            return recompute_compatibility;
        }

        // Recompute the state !
        void recomputeCompatibility() {
            recompute_compatibility = true;
        }

        // Return the state
        int getState() {
            return column_state_;
        }

        // Is column in current solution ?
        bool isInCurrentSolution() {
            return is_in_current_solution_;
        }

        // Set column cost
        void setCost(double cost) {
            cost_ = cost;
        }

        // Get column cost
        double getCost() {
            return cost_;
        }

        // Set column contributions
        void setContribs(std::vector<std::pair<std::string, int>> contribs) {
            contribsIndex_.clear();
            contributions = contribs;
            for (auto pair : contribs) {
                contribsIndex_[pair.first] = pair.second;
            }
        }

        // Find contribution of the column on task taskName
        int findContribution(std::string taskName) {
            if (contribsIndex_.find(taskName) == contribsIndex_.end()) {
                return 0;
            }

            return contribsIndex_[taskName];
        }

        // Change a contribution
        void replaceContrib(std::string oldTask, std::string newTask) {
            if (contribsIndex_.find(oldTask) != contribsIndex_.end()) {
                contribsIndex_[newTask] = contribsIndex_[oldTask];
                contribsIndex_.erase(oldTask);
            }

        }

        // Delete a contributions
        void deleteContrib(std::string taskName) {
            if (contribsIndex_.find(taskName) != contribsIndex_.end()) {
                contribsIndex_.erase(taskName);
            }
        }

        // Get contributions ordered
        std::vector<std::string> getContribsOrder() {
            return contribsOrder;
        }

        // Get contributions non ordered
        std::map<std::string, int> getContribs() {
            return contribsIndex_;
        }

        // Set inc degree of column
        void setPhase(int phase) {
            phase_ = phase;
        }

        // Get inc degree of column
        int getPhase() {
            return phase_;
        }
};
#endif