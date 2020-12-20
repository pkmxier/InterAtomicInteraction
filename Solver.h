#ifndef ATOMICENERGY_SOLVER_H
#define ATOMICENERGY_SOLVER_H
#include "Lattice.h"
#include "Parameters.h"

enum class TableParamType {
    a, E_coh, B, C_11, C_12, C_44, E_sol, E_in, E_on, Size
};

class Solver {
private:
    const double alpha = 1e-6;
    Lattice lattice;
    std::vector<double> target;

    std::pair<int, int> E_in_indices;
    std::pair<Vector3D, Vector3D> E_on_positions;
    double E_coh_A;

    double CalculateLatticeConstant(const Parameters &params, double a_min, double a_max, int n, int itersCount);

    double toEVolt(double energy);
    double Derivative2(double E, const Parameters &params, const std::vector<double> &positiveTransform,
                       const std::vector<double> &negativeTransform);
    std::vector<double> Calculate(const Parameters &params);
    double Eb(int i, const Parameters &params, const std::vector<double> &transform, bool enableCutoff) const;
    double Er(int i, const Parameters &params, const std::vector<double> &transform, bool enableCutoff) const;
public:
    Solver(const Lattice &lattice, const std::vector<double> &target, double energy,
           const std::vector<int> &indices, const std::vector<Vector3D> &positions);

    double MinimizeFunction(Parameters &params);
    double Energy(const Parameters &params, const std::vector<double> &transform = std::vector<double>(), bool enableCutoff = true) const;

    void MoveAtom(int i, const Vector3D &v);
    Lattice & GetLattice();
};


#endif //ATOMICENERGY_SOLVER_H
