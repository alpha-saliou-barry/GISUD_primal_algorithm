#ifndef IB_Column_H
#define IB_Column_H
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <set>

class IB_Column {
    private:
        std::vector<std::string> compatibleBy_;
        std::vector<std::pair<std::string, int>> contributions;
        std::map<std::string, int> contribsIndex_;
        std::set<std::string> tasksInColumn;
        int column_state_;
        bool is_in_current_solution_ = false;
        double cost_;
        int phase_;
        std::string name_;
        bool is_in_p = false;
        bool recompute_compatibility = true;
        std::vector<std::string> contribsOrder;
        bool is_in_p_independent = false;

    public:
        static const int STATE_COMPATIBLE = 1;
        static const int STATE_INCOMPATIBLE = 2;
        static const int STATE_COMPATIBLE_LINEARLY = 3;

        IB_Column (std::string name, std::vector<std::pair<std::string, int>> contribs) : name_(name), contributions(contribs), cost_(0), compatibleBy_() {
            compatibleBy_.clear();
            column_state_ = STATE_INCOMPATIBLE;
            is_in_current_solution_ = false;

            for (auto pair : contribs) {
                contribsIndex_[pair.first] = pair.second;
                tasksInColumn.insert(pair.first);
            }
        };

        bool isSame(std::set<std::string> coveredTasks) {
            return tasksInColumn == coveredTasks;
        }

        void setCompatibleby(std::set<std::string>  cb) {
            compatibleBy_.clear();
            
            for (auto colIndice : cb) {
                compatibleBy_.push_back(colIndice);
            }
        }

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

        std::vector<std::string> getCompatibleBy() {
            return compatibleBy_;
        }


        std::string getName() {
            return name_;
        }
        
        void setName(std::string name) {
            name_ = name;
        }


        void setCompatible() {
            column_state_ = STATE_COMPATIBLE;
        }

        void setIncompatible() {
            column_state_ = STATE_INCOMPATIBLE;
        }

        void setCompatibleLinearly() {
            column_state_ = STATE_COMPATIBLE_LINEARLY;
        }

        void setInCurrentSolution() {
            is_in_current_solution_ = true;
        }

        void setOutCurrentSolution() {
            is_in_current_solution_ = false;
        }

        void changeInP(bool newState) {
            is_in_p = newState;
        }

        bool isInP() {
            return is_in_p;
        }

        void changeInPIndependent(bool new_state) {
            is_in_p_independent = new_state;
        }

        bool isInPIndependent() {
            return is_in_p_independent;
        }

        bool shouldRecomputeCompatibility() {
            return recompute_compatibility;
        }

        void recomputeCompatibility() {
            recompute_compatibility = true;
        }

        int getState() {
            return column_state_;
        }

        bool isInCurrentSolution() {
            return is_in_current_solution_;
        }

        void setCost(double cost) {
            cost_ = cost;
        }

        double getCost() {
            return cost_;
        }

        void setContribs(std::vector<std::pair<std::string, int>> contribs) {
            contribsIndex_.clear();
            contributions = contribs;
            for (auto pair : contribs) {
                contribsIndex_[pair.first] = pair.second;
            }
        }



        int findContribution(std::string taskName) {
            if (contribsIndex_.find(taskName) == contribsIndex_.end()) {
                return 0;
            }

            return contribsIndex_[taskName];
        }

        void replaceContrib(std::string oldTask, std::string newTask) {
            if (contribsIndex_.find(oldTask) != contribsIndex_.end()) {
                contribsIndex_[newTask] = contribsIndex_[oldTask];
                contribsIndex_.erase(oldTask);
            }

        }

        void deleteContrib(std::string taskName) {
            if (contribsIndex_.find(taskName) != contribsIndex_.end()) {
                contribsIndex_.erase(taskName);
            }
        }

        std::vector<std::string> getContribsOrder() {
            return contribsOrder;
        }

        std::map<std::string, int> getContribs() {
            return contribsIndex_;
        }

        void setPhase(int phase) {
            phase_ = phase;
        }

        int getPhase() {
            return phase_;
        }
};
#endif