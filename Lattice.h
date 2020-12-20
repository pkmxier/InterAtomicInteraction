#ifndef ATOMICENERGY_LATTICE_H
#define ATOMICENERGY_LATTICE_H
#include "atom.h"
#include <vector>

class Lattice {
private:
    std::vector<Atom> atoms;
    Vector3D period;
    double multiplier;
public:
    Lattice(const std::vector<Atom> &atoms, Vector3D period, double multiplier);

    int Size() const;
    double Distance(int i, int j, std::vector<double> transform, bool enableCutoff) const;

    void SetPeriod(const Vector3D &p);
    void SetMultiplier(double m);

    Vector3D GetPeriod() const;
    double GetMultiplier() const;
    double Get_r_0() const;

    void PushBack(const Atom &atom);
    void Pop_back();

    static Lattice GenerateLattice(AtomType type, double multiplier);

    Atom & operator [](int i);
    const Atom & operator [](int i) const;
};


#endif //ATOMICENERGY_LATTICE_H
