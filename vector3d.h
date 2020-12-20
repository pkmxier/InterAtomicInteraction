#ifndef VECTOR3D_H
#define VECTOR3D_H
#include <iostream>
#include <vector>

class Vector3D {
public:
    double x, y, z;
    
    Vector3D() {}
    Vector3D(double x, double y, double z);
    Vector3D(const Vector3D &a, const Vector3D &b);
    Vector3D(const std::vector<double> &v);

    double Norm() const;
    friend Vector3D multiply(const Vector3D &lhs, const Vector3D &rhs);         // coordinate-wise multipication

    double & operator[](int i);
    const double & operator[](int i) const;

    Vector3D operator /(double scale);
    friend Vector3D operator +(const Vector3D &lhs, const Vector3D &rhs);
    friend Vector3D operator -(const Vector3D &lhs, const Vector3D &rhs);
    void operator +=(Vector3D rhs);

    friend Vector3D operator *(const Vector3D &lhs, double scale);
    friend Vector3D operator *(double scale, const Vector3D &rhs);

    friend std::ostream & operator <<(std::ostream &os, const Vector3D &v);

    friend bool operator ==(const Vector3D &lhs, const Vector3D &rhs);
    
    static Vector3D i();
    static Vector3D j();
    static Vector3D k();
};



#endif // VECTOR3D_H
