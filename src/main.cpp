// ConsoleApplication1.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include "ISUD.h"
#include <string>
#include <cstring>
#include "CplexMIP.h"
#include <Eigen/Dense>
#include "IB_CompatibilityChecker.h"
#include <numeric>
#include "ProblemReader.h"
#include "InitialSolutions.h"
#include "Files.h"
#include "Params.h"

// Main function
int main(int argc, char* argv[])
{
    std::string path = argv[1];
    char* cplex = argv[2];
    Params params(path + "/params.txt");
    double gap = params.getDoubleParam("gap");

    if(strcmp(cplex, "0") == 0) {
        char* rc = argv[3];
        char* compete = argv[4];
        ISUD_Base problem = constructISUDProblem(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
        ISUD isud(&problem, gap, strcmp(rc, "1") == 0, false, strcmp(compete, "1") == 0);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        isud.solve(path);
        problem.destroyColumns();
    } else {
      ISUD_Base problem = constructISUDProblem(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
      ISUD isud(&problem, gap, false, false, true);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        CplexMIP cplex(&problem);
        std::vector<int> solution;

        double objective = cplex.solve(initialSolution, &solution, false, gap, path);
        problem.destroyColumns();
    }

    return 0;
}
