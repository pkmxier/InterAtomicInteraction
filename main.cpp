#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

#include "Optimizer.h"
#include "Timer.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "Solver.h"
#include <fstream>

#include <omp.h>

namespace pt = boost::property_tree;

void test(const Parameters &params, double energy, double latticeConstant,
          const std::vector<int> &indices, const std::vector<Vector3D> &positions) {
    Vector3D pos1 = Vector3D(0, 0, 0);
    Vector3D pos2 = Vector3D(0.01, 0, 0);
    Vector3D period = Vector3D(1, 1, 1) * std::numeric_limits<double>::max();

    std::ofstream file;
    file.open("../file.py");

    std::vector< std::vector<Atom> > v = {{Atom(AtomType::A, pos1), Atom(AtomType::A, pos2)},
                                          {Atom(AtomType::A, pos1), Atom(AtomType::B, pos2)},
                                          {Atom(AtomType::B, pos1), Atom(AtomType::B, pos2)}};
    std::vector<std::string> names = {"AA", "AB", "BB"};

    Parameters p(std::vector<double>{0.106932, 0.0160873, 1.53027, 10.6992, 2.12453, 2.50854,
                                        0.101786, -0.0288811, 1.60002, 11.9019, 3.84388, 2.73478,
                                        0.0791054, 0.017081, 1.11335, 11.9745, 3.26843, 2.95808,
    });

    std::vector<double> x;

    file << "Energy = [0, 0, 0]\n";
    for (int k = 0; k < v.size(); ++k) {
        Solver solver(Lattice(v[k], period, latticeConstant), std::vector<double>(), energy, indices, positions);

        file << "Energy[" << k << "] = [";
        for (int i = 0; i < 200; ++i) {
            double val = solver.Energy(params, std::vector<double>(), false);
            file << val << ", ";

            if (k == 0)
                x.push_back(solver.GetLattice().Distance(0, 1, std::vector<double>(), false));

            Vector3D newPos = pos2 + Vector3D::i() * i / 100.;
            solver.MoveAtom(1, newPos);
        }
        file << "]\n";
    }

    file << "x = [";
    for (int i = 0; i < x.size(); ++i) {
        file << x[i] << ", ";
    }
    file << "]\n";

    file << "\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "\n"
            "fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey='all')\n"
            "\n"
            "min = 1e9\n"
            "for E in Energy:\n"
            "    cur = np.min(E)\n"
            "    min = cur if cur < min else min\n"
            "\n"
            "names = [\"A-A\", \"A-B\", \"B-B\"]\n"
            "for e, axis, name in zip(Energy, ax, names):\n"
            "    axis.grid(True)\n"
            "    axis.set_ylim(min - 0.3, 5)\n"
            "    axis.set_title(name)\n"
            "    axis.set_xlabel(\"r(Angstrem)\")\n"
            "    axis.set_ylabel(\"Energy(eVolt)\")\n"
            "    axis.axhline(color='red')\n"
            "    axis.axvline(color='red')\n"
            "    axis.plot(x, e)\n"
            "plt.savefig(\"graph.png\")";

    file.close();
}

int main(int argc, char **argv) {
    pt::ptree root;
    pt::read_json("params_V_Ag.json", root);

    std::vector<double> targetParams;
    int n = 18;
    std::pair<Parameters, Parameters> initParams((Parameters(n)), Parameters(n));

    for (auto &it: root.get_child("table")) {
        targetParams.push_back(stod(it.second.data()));
    } // a, E_coh, B, C_11, C_12, C_44, E_sol, E_in, E_on

    int i = 0;
    for (auto &it: root.get_child("init_params")) {
        std::vector<double> values;

        for (auto &it1: root.get_child("init_params." + it.first) ) {
            values.push_back(stod(it1.second.data()));
        }

        for (int j = 0; j < 3; ++j) {
            initParams.first[i + 6 * j] = values[0];
            initParams.second[i + 6 * j] = values[1];
        }

        ++i;
    } // A_0, A_1, ksi, p, q, r_0

    std::vector<int> E_in_indices;
    for (auto &it: root.get_child("E_in_indices")) {
        E_in_indices.push_back(stoi(it.second.data()));
    }

    std::vector<Vector3D> E_on_positions;
    for (auto &it: root.get_child("E_on_positions")) {
        std::vector<double> position;
        for (auto it1: it.second) {
            position.push_back(stod(it1.second.data()));
        }
        E_on_positions.push_back(Vector3D(position));
    }

    double E_coh_A = root.get<double>("E_coh_A");
    double latticeConstant = root.get<double>("latticeConstant");

    Lattice lattice = Lattice::GenerateLattice(AtomType::B, latticeConstant);
    Solver solver(lattice, targetParams, E_coh_A, E_in_indices, E_on_positions);

    Optimizer optimizer(solver, 1, 0.5, 2, 0.5);

    Timer t;
    t.Start();

    Parameters params = optimizer.NelderMeadMinimization(initParams, 1e-6);
    std::cout << "Optimized: \n" << params;

    std::cout << t.ElapsedMilliseconds() << "\n";

    test(params, E_coh_A, latticeConstant, E_in_indices, E_on_positions);
}
