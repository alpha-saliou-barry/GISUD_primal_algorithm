#ifndef IB_CompatibilityChecker_H
#define IB_CompatibilityChecker_H
#include "IB_Column.h"
#include "ISUD_Base.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <set>
#include<Eigen/SparseLU>

class IB_CompatibilityChecker
{
    

private:
    // Pointer on the problem
    ISUD_Base *psolutionMethod_;
    // Final columns candidates
    std::vector<int> finalCandidates;
    // Inverse structure variables
    Eigen::SparseMatrix<float> ipositiveColumnsMatrix_, oposMatrix;
    Eigen::FullPivLU<Eigen::MatrixXf>* inverseStructure_ = NULL, *squareMatrixInverse_ = NULL;

    // Active constraints and positive columns (indices and costs)
    std::vector<std::string> activeConstraints_;
    std::vector<int> indicesActiveConstraints_;
    std::vector<int> subset, posColumns;
    Eigen::VectorXf positiveColumnsCosts_;

    // Other constraints (non active)
    std::vector<int> otherConstraints;

    // Vector of compatibility
    Eigen::VectorXf compatibilityVector;
    std::set<int> positiveColumnsIndices;

    // RandomVector and a matrix
    Eigen::VectorXf randomVector;
    Eigen::MatrixXf theMatrix;

public:
    // Constructor
    IB_CompatibilityChecker(ISUD_Base *psolutionMethod)
    {
        psolutionMethod_ = psolutionMethod;
        squareMatrixInverse_ = NULL;
        inverseStructure_ = NULL;
    }

    // Compute if the row h is dependent
    bool isRowDependent(int i);

    // Compute the linearly independent matrix from positive columns
    void calcIndependentMatrix();

    // Compute the inverse structure
    void calcInverseStructure();

    // Return if column in columns are linear compatible, put the result in solution vector of booleans
    void isLinearlyCompatible(std::vector<IB_Column*>& columns, std::vector<bool>& solution);

    // Get active constraints of the RP
    std::vector<std::string> getActiveConstraints();

    // Get CP active constraints
    std::vector<std::string> getCPActiveConstraints();

    // Convert column to Eigen vector
    Eigen::VectorXf columnToEigenVector(IB_Column* column);

    // Get distance of two points separated by direction of support nonNullColumnsINdices
    int getDistance(std::vector<int> nonNullColumnsIndices);

    // Get incompatible vector and cost for reduction of CP
    void getIncompatibleVectorAndCost(std::vector<IB_Column*> columns, std::vector<Eigen::VectorXf>* vectors, std::vector<double>* colCosts);

    // Set Active constraints
    void setActiveConstraints(std::vector<std::string> ac);

    // Add active constraint in the RP
    void addActiveConstraint(int i) {
        int pos = activeConstraints_.size();
        for (int j = 0; j < activeConstraints_.size(); j++) {
            if (indicesActiveConstraints_[j] > i) {
                pos = j;
                break;
            }
        }

        activeConstraints_.insert(activeConstraints_.begin() + pos, psolutionMethod_->tasks_[i]);
        indicesActiveConstraints_.insert(indicesActiveConstraints_.begin() + pos, i);
    }

    // Get active columns, (positive columns)
    std::set<int> getActiveColumns() {
        std::set<int> ac;
        for (auto i : posColumns) {
            ac.insert(i);
        }

        return ac;
    }

    // Destructor
    ~IB_CompatibilityChecker() {
        if (squareMatrixInverse_ != NULL) {
            delete squareMatrixInverse_;
            squareMatrixInverse_ = NULL;
        }
    }
};
#endif