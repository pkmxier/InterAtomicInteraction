#include "Optimizer.h"
#include <ctime>

Optimizer::Optimizer(const Solver &solver, double a, double b, double g, double s)
    : alpha(a), beta(b), gamma(g), sigma(s), solver(solver) {
    generator = std::mt19937(time(0)); // std::random_device()();
}

Parameters Optimizer::NelderMeadMinimization(std::pair<Parameters, Parameters> &bounds, double epsilon) {
    int n = bounds.first.size();

    std::vector< std::uniform_real_distribution<double> > distributionBounds(n);
    double maxDouble = std::numeric_limits<double>::max();

    for (int i = 0; i < n; ++i) {
        distributionBounds[i] =
                std::uniform_real_distribution<double>(bounds.first[i], std::nextafter(bounds.second[i], maxDouble));
    }

    std::vector<Parameters> points(n + 1, Parameters(n));
    std::vector<double> functionValues(n + 1);

    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < n; ++j) {
            points[i][j] = distributionBounds[j](generator);
        }
    }

    std::pair<int, double> max(-1, std::numeric_limits<double>::lowest());
    std::pair<int, double> mid = max;
    std::pair<int, double> min(-1, std::numeric_limits<double>::max());

    int iteration = 0;
    bool wasShrinked = true;
    double convergence = maxDouble;

    Parameters reflected(n);
    Parameters average(n);
    Parameters expanded(n);
    Parameters shrinked(n);

    do {
        max = std::make_pair(-1, std::numeric_limits<double>::lowest());
        mid = max;

        for (int i = 0; i < n + 1; ++i) {
            if (wasShrinked) {
                functionValues[i] = solver.MinimizeFunction(points[i]);
            }

            if (functionValues[i] > max.second) {
                mid = max;
                max = std::make_pair(i, functionValues[i]);
            } else if (functionValues[i] > mid.second) {
                mid = std::make_pair(i, functionValues[i]);
            } else if (functionValues[i] < min.second) {
                min = std::make_pair(i, functionValues[i]);
            }
        }

        average.FillWith(0);
        for (int i = 0; i < n + 1; ++i) {
            if (i == max.first) {
                continue;
            }

            average += points[i];
        }
        average /= n;

        for (int i = 0; i < n; ++i) {
            reflected[i] = average[i] + alpha * (average[i] - points[max.first][i]);
        }

        double reflectedValue = solver.MinimizeFunction(reflected);
        bool needsCompression = false;
        wasShrinked = false;

        if (reflectedValue < min.second) {
            for (int i = 0; i < n; ++i) {
                expanded[i] = (1 - gamma) * average[i] + gamma * reflected[i];
            }

            double expandedValue = solver.MinimizeFunction(expanded);

            if (expandedValue < reflectedValue) {
                points[max.first] = expanded;
                functionValues[max.first] = expandedValue;
            } else {
                points[max.first] = reflected;
                functionValues[max.first] = reflectedValue;
            }
        } else if (reflectedValue < mid.second) {
            points[max.first] = reflected;
            functionValues[max.first] = reflectedValue;
        } else if (reflectedValue < max.second) {
            points[max.first] = reflected;
            functionValues[max.first] = reflectedValue;
            needsCompression = true;
        } else {
            needsCompression = true;
        }

        if (needsCompression) {
            for (int i = 0; i < n; ++i) {
                shrinked[i] = (1 - beta) * average[i] + beta * points[max.first][i];
            }

            double shrinkedValue = solver.MinimizeFunction(shrinked);

            if (shrinkedValue < max.second) {
                points[max.first] = shrinked;
                functionValues[max.first] = shrinkedValue;
            } else {
                wasShrinked = true;

                for (int i = 0; i < n + 1; ++i) {
                    if (i == min.first) {
                        continue;
                    }

                    for (int j = 0; j < n; ++j) {
                        points[i][j] = points[min.first][j] +
                                       sigma * (points[i][j] - points[min.first][j]);
                    }
                }
            }
        }

        ++iteration;

        convergence = ConvergenceCriterion(points);
        if (iteration % 50 == 0) {
            std::cout << "i = " << iteration << ": err = " << min.second << ", convergence = " << convergence << "\n";
        }
    } while (min.second >= epsilon && convergence >= 1e-3);

    std::cout << "End: i = " << iteration << ": err = " << min.second << ", convergence = " << convergence << "\n";

    return points[min.first];
}

double Optimizer::ConvergenceCriterion(std::vector<Parameters> &points) {
    int n = points.size();
    int m = points[0].size();

    Parameters average(m);
    Parameters dispersion(m);

    for (int i = 0; i < n; ++i) {
        average += points[i];
    }
    average /= n;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            dispersion[j] += (points[i][j] - average[j]) * (points[i][j] - average[j]);
        }
    }

    for (int j = 0; j < m; ++j) {
        dispersion[j] = std::sqrt(dispersion[j] / n);
    }

    std::vector<Parameters> normalizedPoints(points.size(), Parameters(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            normalizedPoints[i][j] = (points[i][j] - average[j]) / dispersion[j];
        }
    }

    double maxNorm = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            maxNorm = std::max(maxNorm, Norm(normalizedPoints[i], normalizedPoints[j]));
        }
    }

    return maxNorm;
}