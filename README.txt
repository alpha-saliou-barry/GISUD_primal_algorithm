GISUD source code

COMPILATION

Code written in C++. Dependencies are CPLEX 20.1 and Eigen 3.4.0 (for linear algebra) : https://eigen.tuxfamily.org/index.php



RUN

In order to run GISUD on an instance located at [INSTANCE_FOLDER] run command :

./main [INSTANCE_FOLDER] [CPLEX] [COLUMN_ADDITION] [DISAGGREGATION] [COMPETITION]

[CPLEX] : 0 if we want to run GISUD and 1 if we want to run CPLEX
[COLUMN_ADDITION] : If GISUD is run, 1 to enable column addition strategy else 0
[DISAGGREGATION] : If GISUD is run, 1 to enable task disaggregation algorithm else 0
[COMPETITION] : If GISUD is run, 1 to enable the comparison between the strategy enabled and ZOOM procedure or 0 otherwise



EXAMPLES :

run only CPLEX on instance "/home/user/instance_test" : 

./main "/home/user/instance_test" 1

run GISUD with column addition strategy (article Primal algorithm with variable addition for the
generalized set-partitioning problem) on instance "/home/user/instance_test" without comparison with ZOOM procedure

./main "/home/user/instance_test" 0 1 0 0