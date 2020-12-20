#ifndef ATOMICENERGY_PARAMETERS_H
#define ATOMICENERGY_PARAMETERS_H
#include <vector>
#include <iostream>

class Parameters {
private:
    std::vector<double> data;
public:
    Parameters(int n);
    Parameters(const std::vector<double> &p);

    int size();
    void FillWith(double value);

    const double & operator [](int i) const;
    double & operator [](int i);

    void operator +=(Parameters &rhs);
    void operator /=(double rhs);

    friend double Norm(const Parameters &a, const Parameters &b);

    friend std::ostream & operator <<(std::ostream &os, const Parameters &p);
};

#endif //ATOMICENERGY_PARAMETERS_H
