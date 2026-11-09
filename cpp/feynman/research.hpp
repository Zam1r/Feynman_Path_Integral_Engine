#pragma once
#include "schrodinger.hpp"
#include <cmath>

namespace feynman {

inline double ho_outside_fixed(double hbar, double wall) {
    if (!(hbar > 0.0) || !(wall > 0.0)) return 0.0;
    return std::erfc(wall / std::sqrt(hbar));
}

inline double ho_outside_turning(double hbar) {
    if (!(hbar > 0.0)) return 0.0;
    return std::erfc(1.0);
}

inline bool double_well_edges(double energy, double a, double b, double& xin, double& xout) {
    if (!(a > 0.0) || !(b > 0.0) || !(energy > 0.0)) return false;
    double barrier = a * b * b * b * b;
    if (!(energy < barrier)) return false;
    double s = std::sqrt(energy / a);
    double inner2 = b * b - s;
    double outer2 = b * b + s;
    if (!(inner2 > 0.0) || !(outer2 > inner2)) return false;
    xin = std::sqrt(inner2);
    xout = std::sqrt(outer2);
    return true;
}

inline double density_between(const Spectrum& s, double lo, double hi) {
    int n = (int)s.x.size();
    if (n < 2 || !(hi > lo)) return 0.0;
    double acc = 0.0;
    for (int i = 0; i + 1 < n; ++i) {
        double xa = s.x[i];
        double xb = s.x[i + 1];
        if (!(xb > xa)) continue;
        double left = xa > lo ? xa : lo;
        double right = xb < hi ? xb : hi;
        if (!(right > left)) continue;
        double ya = s.psi2[i];
        double yb = s.psi2[i + 1];
        double ta = (left - xa) / (xb - xa);
        double tb = (right - xa) / (xb - xa);
        double fa = ya + (yb - ya) * ta;
        double fb = ya + (yb - ya) * tb;
        acc += 0.5 * (fa + fb) * (right - left);
    }
    return acc;
}

}
