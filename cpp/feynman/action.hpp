#pragma once
#include "potential.hpp"
#include <vector>

namespace feynman {

inline double link_action(double left, double right, double eps, double mass, double hbar, const Potential& V) {
    double dx = right - left;
    double mid = 0.5 * (left + right);
    return (mass * dx * dx) / (2.0 * eps * hbar) + (eps / hbar) * V(mid);
}

inline double energy_estimator(const std::vector<double>& x, double eps, double mass, double hbar, const Potential& V) {
    int N = (int)x.size();
    if (N < 1 || !(eps > 0.0) || !(mass > 0.0)) return 0.0;
    double kin = 0.0;
    double pot = 0.0;
    for (int j = 0; j < N; ++j) {
        int k = j + 1;
        if (k == N) k = 0;
        double dx = x[k] - x[j];
        double mid = 0.5 * (x[j] + x[k]);
        kin += (mass * dx * dx) / (2.0 * hbar * eps * eps);
        pot += V(mid) / hbar;
    }
    kin /= (double)N;
    pot /= (double)N;
    return 1.0 / (2.0 * eps) - kin + pot;
}

inline int sign_crossings(const std::vector<double>& x) {
    int N = (int)x.size();
    int flips = 0;
    for (int j = 0; j < N; ++j) {
        int k = j + 1;
        if (k == N) k = 0;
        if (x[j] * x[k] < 0.0) ++flips;
    }
    return flips;
}

inline double periodic_action(const std::vector<double>& x, double eps, double mass, double hbar, const Potential& V) {
    int N = (int)x.size();
    double S = 0.0;
    for (int j = 0; j < N; ++j) {
        int k = j + 1;
        if (k == N) k = 0;
        S += link_action(x[j], x[k], eps, mass, hbar, V);
    }
    return S;
}

}
