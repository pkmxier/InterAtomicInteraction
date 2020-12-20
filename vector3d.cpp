#include "vector3d.h"
#include <cmath>

Vector3D::Vector3D(double x, double y, double z)
    : x(x), y(y), z(z) {
}

Vector3D::Vector3D(const std::vector<double> &v)
    : x(v[0]), y(v[1]), z(v[2]) {
}

double Vector3D::Norm() const {
    return std::sqrt(x * x + y * y + z * z);
}

Vector3D multiply(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}

double operator *(const Vector3D &lhs, const Vector3D &rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vector3D operator *(const Vector3D &lhs, double scale) {
    return Vector3D(lhs.x * scale, lhs.y * scale, lhs.z * scale);
}

Vector3D operator *(double scale, const Vector3D &rhs) {
    return rhs * scale;
}

Vector3D Vector3D::operator /(double scale) {
    return (*this) * (1.0f / scale);
}

void  Vector3D::operator +=(Vector3D rhs) {
    *this = *this + rhs;
}

Vector3D operator ^(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.y * rhs.z - lhs.z * rhs.y,
                    - lhs.x * rhs.z + lhs.z * rhs.x,
                    lhs.x * rhs.y - lhs.y * rhs.x);
}

Vector3D operator +(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

Vector3D operator -(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

std::ostream & operator <<(std::ostream &os, const Vector3D &v) {
    os << v.x << " " << v.y << " " << v.z << std::endl;
    return os;
}

bool operator ==(const Vector3D &lhs, const Vector3D &rhs) {
    double bias = 1e-16;
    return (lhs - rhs).Norm() < bias;
}

Vector3D Vector3D::i() {
    return Vector3D(1, 0, 0);
}

Vector3D Vector3D::j() {
    return Vector3D(0, 1, 0);
}

Vector3D Vector3D::k() {
    return Vector3D(0, 0, 1);
}

double & Vector3D::operator[](int i) {
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            return x;
    }
}

const double & Vector3D::operator[](int i) const {
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            return x;
    }
}
