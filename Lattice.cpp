#include "Lattice.h"
#include <limits>
#include <cmath>

Lattice::Lattice(const std::vector<Atom> &atoms, Vector3D period, double multiplier)
    : atoms(atoms), period(period), multiplier(multiplier) {
}

double Lattice::Distance(int i, int j, std::vector<double> transform, bool enableCutoff) const {
    Vector3D direction = atoms[j].position - atoms[i].position;

    for (int i = 0; i < 3; ++i) {
        if (direction[i] > period[i] / 2) {
            direction[i] -= period[i];
        } else if (direction[i] < -period[i] / 2) {
            direction[i] += period[i];
        }
    }

    if (transform.size() == 3) {
        direction = multiply(direction, Vector3D(transform));
    } else if (transform.size() == 5) {
        direction = Vector3D(direction.x * transform[0] + direction.y * transform[1],
                             direction.x * transform[2] + direction.y * transform[3],
                             direction.z * transform[4]);
    }

    double norm = direction.Norm();
    if (norm > 1.7 && enableCutoff) {
        return std::numeric_limits<double>::max();
    }

    return norm * multiplier;
}

Atom & Lattice::operator[](int i) {
    return atoms[i];
}

const Atom & Lattice::operator[](int i) const {
    return atoms[i];
}

int Lattice::Size() const {
    return atoms.size();
}

void Lattice::SetMultiplier(double m) {
    multiplier = m;
}

void Lattice::SetPeriod(const Vector3D &p) {
    period = p;
}

double Lattice::GetMultiplier() const {
    return multiplier;
}

Vector3D Lattice::GetPeriod() const {
    return period;
}

double Lattice::Get_r_0() const {
    return multiplier / std::sqrt(2);
}

void Lattice::PushBack(const Atom &atom) {
    atoms.push_back(atom);
}

void Lattice::Pop_back() {
    atoms.pop_back();
}

Lattice Lattice::GenerateLattice(AtomType type, double multiplier) {
    std::vector<Atom> basis = {
            Atom(type, Vector3D(0, 0, 0)),
            Atom(type, Vector3D(0.5, 0, 0.5)),
            Atom(type, Vector3D(0.5, 0.5, 0)),
            Atom(type, Vector3D(0, 0.5, 0.5))
    };

    std::vector<Atom> result;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0 ; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (auto atom: basis) {
                    result.push_back(Atom(type, Vector3D(atom.position + Vector3D(i, j, k))));
                }
            }
        }
    }

    return Lattice(result, Vector3D(3, 3, 3), multiplier);
}

