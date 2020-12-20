//
// Created by pkmixer on 17.12.2020.
//

#include "Solver.h"
#include <cmath>
#include <omp.h>

Solver::Solver(const Lattice &lattice, const std::vector<double> &target, double energy,
               const std::vector<int> &indices, const std::vector<Vector3D> &positions)
    : lattice(lattice), target(target), E_coh_A(energy),
      E_in_indices(indices[0], indices[1]),
      E_on_positions(positions[0], positions[1]) {
}

double Solver::Energy(const Parameters &params, const std::vector<double> &transform, bool enableCutoff) const {
    double data[12] = {0,};

    omp_set_dynamic(0);
    #pragma omp parallel num_threads(12)
    {
        int threadsNumber = omp_get_num_threads();
        int thread = omp_get_thread_num();
        int n = lattice.Size();
        int lowerBound = n * thread / double(threadsNumber);
        int upperBound = n * (thread + 1) / double(threadsNumber);

        if (thread == threadsNumber - 1) {
            upperBound = n;
        }

        for (int i = lowerBound; i < upperBound; ++i) {
            data[thread] += Eb(i, params, transform, enableCutoff) + Er(i, params, transform, enableCutoff);
        }
    }

    double result = 0;
    for (int i = 0; i < 12; ++i) {
        result += data[i];
    }
    return result;
}

double Solver::Eb(int i, const Parameters &params, const std::vector<double> &transform, bool enableCutoff) const {
    double result = 0;

    for (int j = 0; j < lattice.Size(); ++j) {
        if (j == i) {
            continue;
        }

        double distance = lattice.Distance(i, j, transform, enableCutoff);
        if (distance > 1e8) {
            continue;
        }

        const double *p = &params[getConnectionTypeIndex(lattice[i], lattice[j])];

        result += p[ParamType::ksi] * p[ParamType::ksi] * std::exp(-2 * p[ParamType::q] * (distance / p[ParamType::r_0] - 1));
    }

    return -std::sqrt(result);
}

double Solver::Er(int i, const Parameters &params, const std::vector<double> &transform, bool enableCutoff) const {
    double result = 0;

    for (int j = 0; j < lattice.Size(); ++j) {
        if (j == i) {
            continue;
        }

        double distance = lattice.Distance(i, j, transform, enableCutoff);
        if (distance > 1e8) {
            continue;
        }

        const double *p = &params[getConnectionTypeIndex(lattice[i], lattice[j])];

        result += (p[ParamType::A_1] * (distance - p[ParamType::r_0]) + p[ParamType::A_0]) * std::exp(-p[ParamType::p] * (distance / p[ParamType::r_0] - 1));
    }

    return result;
}

std::vector<double> Solver::Calculate(const Parameters &params) {
    lattice.SetMultiplier(CalculateLatticeConstant(params, 2, 7, 5, 5));

    double E_0 = Energy(params);
    double E_coh = E_0 / lattice.Size();

    double V_0 = std::pow(lattice.GetMultiplier(), 3) * 3 * 3 * 3;

    double C_11_derivative = Derivative2(E_0, params, {1 + alpha, 1 + alpha, 1}, {1 - alpha, 1 - alpha, 1});
    double C_12_derivative = Derivative2(E_0, params, {1 + alpha, 1 - alpha, 1}, {1 - alpha, 1 + alpha, 1});

    double C_11 = toEVolt((C_11_derivative + C_12_derivative) / (4 * V_0));
    double C_12 = toEVolt((C_11_derivative - C_12_derivative) / (4 * V_0));

    double B_derivative = Derivative2(E_0, params, {1 + alpha, 1 + alpha, 1 + alpha}, {1 - alpha, 1 - alpha, 1 - alpha});
    double B = toEVolt(B_derivative / (9 * V_0));

    double C_44_derivative = Derivative2(E_0, params, {1, alpha, alpha, 1, 1 / (1 - alpha * alpha)},
                                         {1, -alpha, -alpha, 1, 1 / (1 - alpha * alpha)});
    double C_44 = toEVolt(C_44_derivative / (4 * V_0));


    lattice[0].sort = AtomType::A;
    double energyAB = Energy(params);
    double E_sol = energyAB - E_0 - E_coh_A + E_coh;
    lattice[0].sort = AtomType::B;


    Vector3D curPeriod = lattice.GetPeriod();
    lattice.SetPeriod(curPeriod + Vector3D::k() * curPeriod.z);

    double E_surface = Energy(params);

    lattice[E_in_indices.first].sort = AtomType::A;
    double E_adatom_in = Energy(params);

    lattice[E_in_indices.second].sort = AtomType::A;
    double E_dimer_in = Energy(params);

    double E_in = E_dimer_in - E_surface - 2 * (E_adatom_in - E_surface);
    lattice[E_in_indices.first].sort = AtomType::B;
    lattice[E_in_indices.second].sort = AtomType::B;


    lattice.PushBack(Atom(AtomType::A, E_on_positions.first));
    double E_adatom_on = Energy(params);

    lattice.PushBack(Atom(AtomType::A, E_on_positions.second));
    double E_dimer_on = Energy(params);

    double E_on = E_dimer_on - E_surface - 2 * (E_adatom_on - E_surface);

    lattice.Pop_back();
    lattice.Pop_back();

    lattice.SetPeriod(curPeriod);

    std::vector<double> tableParams(static_cast<int>(TableParamType::Size));

    tableParams[static_cast<int>(TableParamType::a)] = lattice.GetMultiplier(); //TODO: calculate
    tableParams[static_cast<int>(TableParamType::E_coh)] = E_coh;
    tableParams[static_cast<int>(TableParamType::C_11)] = C_11;
    tableParams[static_cast<int>(TableParamType::C_12)] = C_12;
    tableParams[static_cast<int>(TableParamType::B)] = B;
    tableParams[static_cast<int>(TableParamType::C_44)] = C_44;
    tableParams[static_cast<int>(TableParamType::E_sol)] = E_sol;
    tableParams[static_cast<int>(TableParamType::E_in)] = E_in;
    tableParams[static_cast<int>(TableParamType::E_on)] = E_on;

    return tableParams;
}

double Solver::toEVolt(double energy) {
    return 1.602 * energy;
}

double Solver::MinimizeFunction(Parameters &params) {
    double result = 0;

    std::vector<double> predicted = Calculate(params);

    for (int i = 0; i < predicted.size(); ++i) {
        double error = predicted[i] / target[i] - 1;

        result += error * error;
    }

    return result / (predicted.size() - 1);
}

void Solver::MoveAtom(int i, const Vector3D &v) {
    lattice[i].position = v;
}

double Solver::Derivative2(double E, const Parameters &params, const std::vector<double> &positiveTransform,
                           const std::vector<double> &negativeTransform) {
    double positiveEnergy = Energy(params, positiveTransform);
    double negativeEnergy = Energy(params, negativeTransform);

    return (positiveEnergy - 2 * E + negativeEnergy) / (alpha * alpha);
}

Lattice &Solver::GetLattice() {
    return lattice;
}

double Solver::CalculateLatticeConstant(const Parameters &params, double a_min, double a_max, int n, int itersCount) {
    std::pair<double, double> min(0, 1e30);
    double width;

    int iter = 0;
    do {
        width = a_max - a_min;
        min = std::make_pair(0, 1e8);

        for (int i = 0; i < n; ++i) {
            double point = a_min + i * width / n;

            lattice.SetMultiplier(point);
            double value = Energy(params);

            if (value < min.second) {
                min = std::make_pair(point, value);
            }
        }

        width /= n;
        a_min = min.first - width / 2;
        a_max = min.first + width / 2;

        ++iter;
    } while(iter < itersCount);

    return min.first;
}
