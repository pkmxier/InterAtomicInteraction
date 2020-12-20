#include "atom.h"

Atom::Atom(AtomType sort, const Vector3D &pos)
    : sort(sort), position(pos) {
}

int getConnectionTypeIndex(Atom lhs, Atom rhs) {
    AtomType a = lhs.sort;
    AtomType b = rhs.sort;

    ConnectionType type;

    if (a == AtomType::A && b == AtomType::B || a == AtomType::B && b == AtomType::A) {
        type = ConnectionType::AB;
    } else if (a == AtomType::A && b == AtomType::A) {
        type = ConnectionType::AA;
    } else {
        type = ConnectionType::BB;
    }

    return ParamType::Size * type;
}