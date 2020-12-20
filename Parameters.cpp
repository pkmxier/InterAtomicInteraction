#include "Parameters.h"
#include <cmath>
#include "atom.h"

Parameters::Parameters(int n)
    : data(n) {
}

Parameters::Parameters(const std::vector<double> &p)
    : data(p) {
}

int Parameters::size() {
    return data.size();
}

void Parameters::FillWith(double value) {
    std::fill(data.begin(), data.end(), value);
}

const double & Parameters::operator [](int i) const {
    return data[i];
}

double & Parameters::operator [](int i) {
    return data[i];
}

void Parameters::operator +=(Parameters &rhs) {
    int n = data.size();
    for (int i = 0; i < n; ++i) {
        data[i] += rhs.data[i];
    }
}

void Parameters::operator /=(double rhs) {
    int n = data.size();
    for (int i = 0; i < n; ++i) {
        data[i] /= rhs;
    }
}

double Norm(const Parameters &a, const Parameters &b) {
    double norm = 0;

    for (int i = 0; i < a.data.size(); ++i) {
        norm += (a[i] - b[i]) * (a[i] - b[i]);
    }

    return std::sqrt(norm);
}

std::ostream & operator <<(std::ostream &os, const Parameters &p) {
    os << "{";
    for (int i = 0; i < p.data.size(); ++i) {
        os << p.data[i] << ", ";

        if ((i + 1) % ParamType::Size == 0) {
            os << std::endl;
        }
    }
    os << "}" << std::endl;

    return os;
}