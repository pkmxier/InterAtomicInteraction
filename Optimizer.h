#ifndef ATOMICENERGY_OPTIMIZER_H
#define ATOMICENERGY_OPTIMIZER_H
#include <random>
#include "Solver.h"

class Optimizer {
private:
    double alpha, beta, gamma, sigma;
    std::mt19937 generator;
    Solver solver;

public:
    Optimizer(const Solver &solver, double alpha = 1, double beta = 0.5, double gamma = 1, double sigma = 0.5);

    Parameters NelderMeadMinimization(std::pair<Parameters, Parameters> &bounds, double epsilon);
    double ConvergenceCriterion(std::vector<Parameters> &points);
};


#endif //ATOMICENERGY_OPTIMIZER_H
