#ifndef __ATOM_H__
#define __ATOM_H__
#include "vector3d.h"

enum AtomType {
    A, B
};

enum ConnectionType {
    AA, AB, BB
};

enum ParamType {
    A_0, A_1, ksi, p, q, r_0, Size
};

struct Atom {
    AtomType sort;
    Vector3D position;
    
    Atom(AtomType sort, const Vector3D &pos);

    friend int getConnectionTypeIndex(Atom lhs, Atom rhs);
};

#endif //__ATOM_H__
