// ConsoleApplication1.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include "ISUD.h"
#include <string>
#include <ilcplex/ilocplex.h>
#include <cstring>
#include "CplexMIP.h"
#include <Eigen/Dense>
#include "IB_CompatibilityChecker.h"
#include <numeric>
#include "Problems.h"

// Main function
int main(int argc, char* argv[])
{
  //modifyRhs();
  //return 0;
  //modifyRhs("/home/barralph/Documents/TestsArticles/vcsMB50_1", "/home/barralph/Documents/TestsArticles/vcsMB50_2");
  //generateProblems();
    //computeCplexSolutions();
  //return 0;
    std::string path = argv[1];
    char* cplex = argv[2];
    if(strcmp(cplex, "0") == 0) {
        char* rc = argv[3];
        //char* dis = argv[4];
        char* compete = argv[4];
        ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");
        ISUD isud(&problem, strcmp(rc, "1") == 0, false, strcmp(compete, "1") == 0);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        isud.solve(path);
        problem.destroyColumns();
    } else {
      ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt","",
                        path + "/fixed_cost.txt");
        ISUD isud(&problem, false, false, true);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        CplexMIP cplex(&problem);
        std::vector<int> solution;

        double objective = cplex.solve(initialSolution, &solution, false, path);
        problem.destroyColumns();
    }

    return 0;
}
