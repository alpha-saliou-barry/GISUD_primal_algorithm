#include "IB_CompatibilityChecker.h"


typedef Eigen::Triplet<double> T;

// Convert column to Eigen vector
Eigen::VectorXf IB_CompatibilityChecker::columnToEigenVector(IB_Column* column) {
    Eigen::VectorXf vector(psolutionMethod_->tasks_.size());
    for(int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        vector(i) = 0;
    }


    std::map<std::string, int> contribs;
    for(auto contrib: column->getContribs()) {
        std::string task = contrib.first;
        int contrib_value = contrib.second;

        contribs[task] = contrib_value;
    }

    for(int i = 0; i < psolutionMethod_->tasks_.size();i++) {
        if(contribs.find(psolutionMethod_->tasks_[i]) != contribs.end()) {
            vector(i) = contribs[psolutionMethod_->tasks_[i]];
        }
    }


    return vector;
}

// Set Active constraints
void IB_CompatibilityChecker::setActiveConstraints(std::vector<std::string> ac) {
    activeConstraints_ = ac;
    indicesActiveConstraints_.clear();
    for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        if (std::find(activeConstraints_.begin(), activeConstraints_.end(), psolutionMethod_->tasks_[i]) != activeConstraints_.end()) {
            indicesActiveConstraints_.push_back(i);
        }
    }
}

// Get distance of two points separated by direction of support nonNullColumnsINdices
int IB_CompatibilityChecker::getDistance(std::vector<int> nonNullColumnsIndices) {
    Eigen::MatrixXf matrix(psolutionMethod_->tasks_.size(), nonNullColumnsIndices.size());
    for (int i = 0; i < nonNullColumnsIndices.size(); i++) {
        matrix(Eigen::seq(0, Eigen::last), i) = columnToEigenVector(psolutionMethod_->columns_[nonNullColumnsIndices[i]]);
    }

    return nonNullColumnsIndices.size() - matrix.colPivHouseholderQr().rank();
}

// Compute if the row h is dependent
bool IB_CompatibilityChecker::isRowDependent(int h) {
    Eigen::VectorXf rhs = Eigen::Map<Eigen::VectorXi>(psolutionMethod_->rhs_.data(), psolutionMethod_->rhs_.size()).cast<float>();
    std::cout << theMatrix.transpose()(Eigen::seq(0, Eigen::last), indicesActiveConstraints_).fullPivLu().rank() << " | " << indicesActiveConstraints_.size() << std::endl;
    Eigen::VectorXf coef = theMatrix.transpose()(Eigen::seq(0, Eigen::last), indicesActiveConstraints_).fullPivLu().solve(theMatrix(h, Eigen::seq(0, Eigen::last)).transpose());
    for (int i = 0; i < coef.size(); i++) {
        if (coef(i) > 0) {
            int constraint = indicesActiveConstraints_[i];
            int sum = 0;
            for (int j = 0; j < theMatrix.cols(); j++) {
                if (theMatrix(constraint, j) == 1 && psolutionMethod_->columns_[finalCandidates[j]]->isInCurrentSolution()) {
                    std::cout << "J: " << finalCandidates[j] << std::endl;
                    sum += 1;
                }
            }

            std::cout << "The sum is : " << sum << " as rhs " << rhs[i] << std::endl;
        }
    }
    std::cout << theMatrix(h, Eigen::seq(0, Eigen::last)) << std::endl;
    std::cout << "Rhs : " << rhs(h) << " | " << rhs(indicesActiveConstraints_).dot(coef) << std::endl;
    
    return true;
}

// Compute the linearly independent matrix from positive columns
void IB_CompatibilityChecker::calcIndependentMatrix() {
    positiveColumnsIndices.clear();
    std::vector<IB_Column*> positiveColumns;
    std::vector<int> indicesPosColumns;
    for(int i = 0; i < psolutionMethod_->columns_.size();i++) {
        IB_Column* column = psolutionMethod_->columns_[i];
        if(column->isInP()) {
            positiveColumns.push_back(column);
            indicesPosColumns.push_back(i);
        }
    }


    std::map<std::string, int> tasksIndices;
    for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        tasksIndices[psolutionMethod_->tasks_[i]] = i;
    }
    
    // Extraction d'un ensemble de colonnes linéairement indépendant maximal
    int matrix_size = 0;
    Eigen::SparseMatrix<float> current_matrix(psolutionMethod_->tasks_.size(), positiveColumns.size());
    std::vector<float> costs;
    std::vector<T> triplets;

    std::cout << positiveColumns.size() << " pos cols" << std::endl;
    for(int j = 0; j < positiveColumns.size(); j++) {
        for (auto contrib : positiveColumns[j]->getContribs()) {
            if (contrib.second != 0) {
                triplets.push_back(T(tasksIndices[contrib.first], j, contrib.second));
            }
        }
    }

    current_matrix.setFromTriplets(triplets.begin(), triplets.end());
    //Décomposition LU
    Eigen::VectorXi independentRows2(psolutionMethod_->tasks_.size()), independentColumns2(positiveColumns.size());
    for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        independentRows2(i) = 0;
    }

    for (int i = 0; i < positiveColumns.size(); i++) {
        independentColumns2(i) = 0;
    }


    Eigen::FullPivLU<Eigen::MatrixXf> lu = current_matrix.toDense().fullPivLu();
    Eigen::VectorXf diagonal = lu.matrixLU().triangularView<Eigen::Upper>().toDenseMatrix().diagonal();
    for (int i = 0; i < diagonal.size(); i++) {
        if (fabs(diagonal[i]) > 1e-3) {
            independentRows2[i] = 1;
            independentColumns2[i] = 1;
        }
        
    }

    Eigen::VectorXi independentRowsEigen = lu.permutationP().inverse() * independentRows2;
    Eigen::VectorXi independentColumnsEigen = lu.permutationQ() * independentColumns2;

    activeConstraints_.clear();
    indicesActiveConstraints_.clear();
    otherConstraints.clear();
    for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        if (independentRowsEigen[i] > 1e-4) {
            activeConstraints_.push_back(psolutionMethod_->tasks_[i]);
            indicesActiveConstraints_.push_back(i);
        }
        else {
            otherConstraints.push_back(i);
        }
    }

    subset.clear();
    posColumns.clear();
    for (int i = 0; i < positiveColumns.size(); i++) {
        if (independentColumnsEigen[i] > 1e-4) {
            positiveColumns[i]->changeInPIndependent(true);
            posColumns.push_back(indicesPosColumns[i]);
            subset.push_back(i);
            positiveColumnsIndices.insert(indicesPosColumns[i]);
        }
        else {
            positiveColumns[i]->changeInPIndependent(false);
        }
    }


    positiveColumnsCosts_ = Eigen::VectorXf(subset.size());
    for (int i = 0; i < subset.size();i++) {
        positiveColumnsCosts_(i) = positiveColumns[subset[i]]->getCost();
    }

    ipositiveColumnsMatrix_ = current_matrix.toDense()(Eigen::seq(0, Eigen::last), subset).sparseView();
}

// Compute the inverse structure
void IB_CompatibilityChecker::calcInverseStructure() {
    Eigen::SparseMatrix<float> current_matrix = ipositiveColumnsMatrix_.toDense()(indicesActiveConstraints_, Eigen::seq(0, Eigen::last)).sparseView();
    std::cout << activeConstraints_.size() << " active constraints." << std::endl;
    std::cout << current_matrix.cols() << " positive independent columns." << std::endl;
    squareMatrixInverse_ = new Eigen::FullPivLU<Eigen::MatrixXf>();
    squareMatrixInverse_->compute(current_matrix);
    //std:: << squareMatrixInverse_->rank() << " rang." << std::endl;
    Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>> transposeSquareMatrixInverse(current_matrix.transpose());

    std::cout << "Rank : " << transposeSquareMatrixInverse.rank() << std::endl;
    std::cout << "Current matrix inversed." << std::endl;
    
    oposMatrix = Eigen::SparseMatrix<float>(otherConstraints.size(), ipositiveColumnsMatrix_.cols());
    std::vector<T> triplets;
    for (int i = 0; i < otherConstraints.size(); i++) {
        for (int j = 0; j < ipositiveColumnsMatrix_.cols(); j++) {
            if (ipositiveColumnsMatrix_.coeff(otherConstraints[i], j) > 1e-4) {
                triplets.push_back(T(i, j, ipositiveColumnsMatrix_.coeff(otherConstraints[i], j)));
            }
        }
    }

    oposMatrix.setFromTriplets(triplets.begin(), triplets.end());
    randomVector = Eigen::VectorXf::Random(oposMatrix.rows());
    compatibilityVector = transposeSquareMatrixInverse.solve(oposMatrix.transpose() * randomVector);
}

// Get active constraints of the RP
std::vector<std::string> IB_CompatibilityChecker::getActiveConstraints() {
    return activeConstraints_;
}

// Return if column in columns are linear compatible, put the result in solution vector of booleans
void  IB_CompatibilityChecker::isLinearlyCompatible(std::vector<IB_Column*>& columns, std::vector<bool>& solution) {
    solution.clear();
    std::map<std::string, int> tasksCorrespondance;
    for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        tasksCorrespondance[psolutionMethod_->tasks_[i]] = i;
    }

    Eigen::SparseMatrix<float> matrix(psolutionMethod_->tasks_.size(), columns.size());
    Eigen::SparseMatrix<float> candidateMatrix1(indicesActiveConstraints_.size(), columns.size()), candidateMatrix2(otherConstraints.size(), columns.size());
    std::map<int, int> constraintsIndexes, constraintsIndexes2;

    for (int i = 0; i < indicesActiveConstraints_.size(); i++) {
        constraintsIndexes[indicesActiveConstraints_[i]] = i;
    }

    for (int i = 0; i < otherConstraints.size(); i++) {
        constraintsIndexes2[otherConstraints[i]] = i;
    }

    std::vector<T> triplets, triplets2, triplets3;
    for (int i = 0; i < columns.size();i++) {
        IB_Column* column = columns[i];
        for (auto contrib : column->getContribs()) {
            if (contrib.second != 0) {
                if (constraintsIndexes.find(tasksCorrespondance[contrib.first]) != constraintsIndexes.end()) {
                    triplets2.push_back(T(constraintsIndexes[tasksCorrespondance[contrib.first]], i, contrib.second));
                }
                else {
                    triplets3.push_back(T(constraintsIndexes2[tasksCorrespondance[contrib.first]], i, contrib.second));
                }

                triplets.push_back(T(tasksCorrespondance[contrib.first], i, contrib.second));
            }
        }
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());
    candidateMatrix1.setFromTriplets(triplets2.begin(), triplets2.end());
    candidateMatrix2.setFromTriplets(triplets3.begin(), triplets3.end());

    /*
    
    Eigen::MatrixXf solution2 = inverseStructure_->solve(matrix.toDense());
    std::cout << "solution calculee" << std::endl;
    Eigen::MatrixXf product = ipositiveColumnsMatrix_ * solution2.sparseView();
    Eigen::MatrixXf matrix2 = matrix.toDense();
    for (int i = 0; i < columns.size(); i++) {
        solution.push_back(matrix2(Eigen::seq(0, Eigen::last), i).isApprox(product(Eigen::seq(0, Eigen::last), i), 1e-4));
    }*/
    
    

    Eigen::VectorXf solution2 = candidateMatrix1.transpose() * compatibilityVector - candidateMatrix2.transpose() * randomVector;
    std::vector<int> candidateCols;
    for (int i = 0; i < columns.size(); i++) {
        if (fabs(solution2(i)) < 1e-2) {
            
            solution.push_back(true);
            candidateCols.push_back(i);
        }
        else {
            if (columns[i]->isInCurrentSolution()) {
                std::cout << solution2(i) << " in the current solution" << std::endl;
            }

            solution.push_back(false);
        }
    }
    
    /*
    std::vector<int> candidateCols;

    for (int i = 0; i < columns.size(); i++) {
        solution.push_back(true);
        candidateCols.push_back(i);
    }*/
    std::cout << candidateCols.size() << " candidate columns." << std::endl;

    int n = 0;
    Eigen::SparseMatrix<float> candidateMatrix(indicesActiveConstraints_.size(), candidateCols.size()), candidateMatrixO(otherConstraints.size(),
        candidateCols.size());
    std::vector<T> tripletsc, tripletscO;

    for (int i = 0; i < candidateCols.size(); i++) {
        IB_Column* column = columns[candidateCols[i]];
        for (auto contrib : column->getContribs()) {
            if (contrib.second != 0) {
                if (constraintsIndexes.find(tasksCorrespondance[contrib.first]) != constraintsIndexes.end()) {
                    tripletsc.push_back(T(constraintsIndexes[tasksCorrespondance[contrib.first]], i, contrib.second));
                }
                else {
                    tripletscO.push_back(T(constraintsIndexes2[tasksCorrespondance[contrib.first]], i, contrib.second));
                }
            }
        }
    }

    candidateMatrix.setFromTriplets(tripletsc.begin(), tripletsc.end());
    candidateMatrixO.setFromTriplets(tripletscO.begin(), tripletscO.end());


    Eigen::MatrixXf solutionMatrix = oposMatrix * squareMatrixInverse_->solve(candidateMatrix.toDense()) - candidateMatrixO;
    
    finalCandidates.clear();

    std::set<int> pcolsTreated;

    for (int i = 0; i < candidateCols.size(); i++) {
        if (fabs(solutionMatrix(Eigen::seq(0, Eigen::last), i).maxCoeff()) > 1e-4 || fabs(solutionMatrix(Eigen::seq(0, Eigen::last), i).minCoeff()) > 1e-4) {
            if (columns[candidateCols[i]]->isInP()) {
                std::cout << "One P column non compatible" << std::endl;
                pcolsTreated.insert(candidateCols[i]);
                std::cout << "Coefficients : " << fabs(solutionMatrix(Eigen::seq(0, Eigen::last), i).maxCoeff()) << " " << fabs(solutionMatrix(Eigen::seq(0, Eigen::last), i).minCoeff()) << std::endl;
            }

            solution[candidateCols[i]] = false;
            n += 1;
        }
        else {
            finalCandidates.push_back(candidateCols[i]);
        }
    }

    std::cout << n << " delete columns." << std::endl;
}

// Get incompatible vector and cost for reduction of CP
void IB_CompatibilityChecker::getIncompatibleVectorAndCost(std::vector<IB_Column*> columns, std::vector<Eigen::VectorXf>* vectors, std::vector<double>* colCosts) {
    std::map<std::string, int> tasksCorrespondance;
    for (int i = 0; i < psolutionMethod_->tasks_.size(); i++) {
        tasksCorrespondance[psolutionMethod_->tasks_[i]] = i;
    }

    Eigen::SparseMatrix<float> matrix(psolutionMethod_->tasks_.size(), columns.size());
    std::vector<T> triplets;
    Eigen::VectorXf pcolCosts = Eigen::Map<Eigen::VectorXf>(positiveColumnsCosts_.data(), positiveColumnsCosts_.size());
    Eigen::VectorXf pastColCosts(columns.size());
    for (int i = 0; i < columns.size(); i++) {
        IB_Column* column = columns[i];
        for (auto contrib : column->getContribs()) {
            if (contrib.second != 0) {
                triplets.push_back(T(tasksCorrespondance[contrib.first], i, contrib.second));
            }
        }

        pastColCosts(i) = (column->getCost());
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::MatrixXf denseMatrix = matrix.toDense();
    
    

    //Eigen::MatrixXf coefficients = squareMatrixInverse_->solve(matrix.toDense()(indicesActiveConstraints_, Eigen::seq(0, Eigen::last)));
    Eigen::SparseMatrix<float> coefficients2 = (squareMatrixInverse_->permutationP() * denseMatrix(indicesActiveConstraints_, Eigen::seq(0, Eigen::last))).sparseView();
    Eigen::SparseMatrix<float> lu_matrix = squareMatrixInverse_->matrixLU().sparseView();
    Eigen::SparseMatrix<float> sparseL = (squareMatrixInverse_->matrixLU().triangularView<Eigen::StrictlyLower>().toDenseMatrix() + Eigen::MatrixXf::Identity(activeConstraints_.size(), activeConstraints_.size())).sparseView();
    const Eigen::TriangularView<const Eigen::SparseMatrix<float>, Eigen::Lower> sparseLTriangular = sparseL.triangularView<Eigen::Lower>();
    sparseLTriangular.solveInPlace(coefficients2);
    Eigen::MatrixXf coefficients = squareMatrixInverse_->permutationQ() * lu_matrix.triangularView<Eigen::Upper>().solve(coefficients2.toDense());
    Eigen::MatrixXf vectorsEigen = denseMatrix(otherConstraints, Eigen::seq(0, Eigen::last)) - oposMatrix * coefficients;
    

    Eigen::VectorXf colCostsEigen = pastColCosts - coefficients.transpose() * pcolCosts;

    for (int i = 0; i < vectorsEigen.cols(); i++) {
        vectors->push_back(vectorsEigen(Eigen::seq(0, Eigen::last), i));
        colCosts->push_back(colCostsEigen[i]);
    }
}

// Get CP active constraints
std::vector<std::string> IB_CompatibilityChecker::getCPActiveConstraints() {
    std::vector<std::string> activeConstraints;
    for (int i = 0; i < otherConstraints.size(); i++) {
        activeConstraints.push_back(psolutionMethod_->tasks_[otherConstraints[i]]);
    }

    return activeConstraints;
}