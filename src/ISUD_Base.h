#ifndef ISUD_BASE_H
#define ISUD_BASE_H
#include <vector>
#include "IB_Column.h"


class ISUD_Base {
    public:
        std::vector< IB_Column* > columns_;
        std::vector< std::string > tasks_;
        std::vector <int> rhs_;
        std::map<std::string, int> colsIndices;
        std::map<std::string, std::vector<std::string>> columnsPerTask;
        double fixed_cost_;
        int sum_bi = 0;

        ISUD_Base(std::vector<std::string> tasks, std::vector <int> rhs, std::vector<IB_Column*> columns, double fixed_cost = 0) : columns_(columns), tasks_(tasks), rhs_(rhs),
        fixed_cost_(fixed_cost) {
            for (int i = 0; i < columns.size(); i++) {
                colsIndices[columns[i]->getName()] = i;
                for (auto contrib : columns[i]->getContribs()) {
                    if (!columnsPerTask.count(contrib.first)) {
                        columnsPerTask[contrib.first] = std::vector<std::string>();
                    }

                    columnsPerTask[contrib.first].push_back(columns[i]->getName());
                }
            }

            for(int i = 0; i < rhs.size();i++) {
                sum_bi += rhs[i];
            }
        }

        bool hasColumn(std::string columnName) {
            return colsIndices.find(columnName) != colsIndices.end();
        }

        void addColumn(IB_Column* column) {
            columns_.push_back(column);
            colsIndices[column->getName()] = columns_.size() - 1;
            for (auto contrib : column->getContribs()) {
                if (!columnsPerTask.count(contrib.first)) {
                    columnsPerTask[contrib.first] = std::vector<std::string>();
                }

                columnsPerTask[contrib.first].push_back(column->getName());
            }
        }

        int addArtificialColumn(IB_Column* column) {
            columns_.push_back(column);
            colsIndices[column->getName()] = columns_.size()-1;
            return columns_.size() - 1;
        }

        void removeColumn(int column_id) {
            colsIndices.erase(columns_[column_id]->getName());
            columns_.erase(columns_.begin() + column_id);
        }

        int getColumnIndice(std::string name) {
            return colsIndices[name];
        }

        void destroyColumns() {
            for (int i = 0; i < columns_.size(); i++) {
                delete columns_[i];
                columns_[i] = NULL;
            }
        }
};

#endif