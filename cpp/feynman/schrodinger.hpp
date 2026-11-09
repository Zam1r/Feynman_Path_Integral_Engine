#pragma once
#include "potential.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

namespace feynman {

struct Spectrum {
    std::vector<double> x;
    std::vector<double> psi;
    std::vector<double> psi2;
    double energy;
    double residual;
    bool ok;
};

inline int sturm_count(const std::vector<double>& d, const std::vector<double>& off, double lambda) {
    int n = (int)d.size();
    int count = 0;
    double u = d[0] - lambda;
    if (u < 0.0) ++count;
    for (int i = 1; i < n; ++i) {
        if (u > -1e-30 && u < 1e-30) u = (u < 0.0) ? -1e-30 : 1e-30;
        u = (d[i] - lambda) - (off[i - 1] * off[i - 1]) / u;
        if (u < 0.0) ++count;
    }
    return count;
}

inline bool thomas(const std::vector<double>& diag, const std::vector<double>& off, const std::vector<double>& rhs, std::vector<double>& y) {
    int n = (int)diag.size();
    std::vector<double> cp(n), rp(n);
    if (!(std::abs(diag[0]) > 1e-18)) return false;
    double super0 = (n > 1) ? off[0] : 0.0;
    cp[0] = super0 / diag[0];
    rp[0] = rhs[0] / diag[0];
    for (int i = 1; i < n; ++i) {
        double sup = (i < n - 1) ? off[i] : 0.0;
        double sub = off[i - 1];
        double den = diag[i] - sub * cp[i - 1];
        if (!(std::abs(den) > 1e-18)) return false;
        cp[i] = sup / den;
        rp[i] = (rhs[i] - sub * rp[i - 1]) / den;
    }
    y[n - 1] = rp[n - 1];
    for (int i = n - 2; i >= 0; --i) y[i] = rp[i] - cp[i] * y[i + 1];
    return true;
}

inline Spectrum solve_schrodinger(const Potential& V, double mass, double hbar, double xmin, double xmax, int n, int level = 0) {
    Spectrum spec;
    spec.energy = 0.0;
    spec.residual = 0.0;
    spec.ok = false;
    if (n < 8 || !(xmax > xmin) || !(mass > 0.0) || !(hbar > 0.0)) return spec;
    double dx = (xmax - xmin) / (double)(n + 1);
    double diag_kin = (hbar * hbar) / (mass * dx * dx);
    double off_val = -0.5 * diag_kin;
    std::vector<double> d(n), off(n - 1, off_val), x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = xmin + (double)(i + 1) * dx;
        d[i] = diag_kin + V(x[i]);
    }
    double vmin = d[0] - diag_kin;
    for (int i = 1; i < n; ++i) vmin = std::min(vmin, d[i] - diag_kin);
    int need = level + 1;
    if (need < 1) need = 1;
    double lo = vmin - 2.0;
    double hi = vmin + 0.25;
    for (int guard = 0; guard < 80 && sturm_count(d, off, hi) < need; ++guard) {
        hi = vmin + std::ldexp(1.0, guard);
    }
    for (int guard = 0; guard < 80 && sturm_count(d, off, lo) >= need; ++guard) lo -= 2.0;
    for (int it = 0; it < 80; ++it) {
        if (hi - lo < 1e-12 * (1.0 + std::abs(hi))) break;
        double mid = 0.5 * (lo + hi);
        if (sturm_count(d, off, mid) >= need) hi = mid;
        else lo = mid;
    }
    double lambda = lo - 1e-8 * (1.0 + std::abs(lo));
    std::vector<double> diag(n), rhs(n, 1.0), y(n, 0.0);
    bool solved = false;
    for (int nudge = 0; nudge < 6 && !solved; ++nudge) {
        double shift = lambda - (double)nudge * 1e-8 * (1.0 + std::abs(lambda));
        for (int i = 0; i < n; ++i) diag[i] = d[i] - shift;
        for (int i = 0; i < n; ++i) {
            rhs[i] = (level == 0) ? 1.0 : std::sin((double)(level + 1) * 3.141592653589793 * (double)(i + 1) / (double)(n + 1));
        }
        bool step_ok = true;
        for (int iter = 0; iter < 40; ++iter) {
            if (!thomas(diag, off, rhs, y)) {
                step_ok = false;
                break;
            }
            double norm2 = 0.0;
            for (int i = 0; i < n; ++i) norm2 += y[i] * y[i];
            if (!(norm2 > 0.0) || !std::isfinite(norm2)) {
                step_ok = false;
                break;
            }
            double inv = 1.0 / std::sqrt(norm2);
            for (int i = 0; i < n; ++i) rhs[i] = y[i] * inv;
        }
        solved = step_ok;
    }
    if (!solved) return spec;
    double ray_num = 0.0;
    double ray_den = 0.0;
    double res2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double Hv = d[i] * rhs[i];
        if (i > 0) Hv += off[i - 1] * rhs[i - 1];
        if (i + 1 < n) Hv += off[i] * rhs[i + 1];
        ray_num += rhs[i] * Hv;
        ray_den += rhs[i] * rhs[i];
    }
    double energy = ray_num / ray_den;
    for (int i = 0; i < n; ++i) {
        double Hv = d[i] * rhs[i];
        if (i > 0) Hv += off[i - 1] * rhs[i - 1];
        if (i + 1 < n) Hv += off[i] * rhs[i + 1];
        double r = Hv - energy * rhs[i];
        res2 += r * r;
    }
    double mass_psi = 0.0;
    for (int i = 0; i < n; ++i) mass_psi += rhs[i] * rhs[i];
    mass_psi *= dx;
    if (!(mass_psi > 0.0)) return spec;
    double scale = 1.0 / std::sqrt(mass_psi);
    double sign = 0.0;
    for (int i = 0; i < n; ++i) sign += rhs[i];
    if (sign < 0.0) scale = -scale;
    spec.x.resize(n);
    spec.psi.resize(n);
    spec.psi2.resize(n);
    for (int i = 0; i < n; ++i) {
        spec.x[i] = x[i];
        spec.psi[i] = rhs[i] * scale;
        spec.psi2[i] = spec.psi[i] * spec.psi[i];
    }
    spec.energy = energy;
    spec.residual = std::sqrt(res2 / ray_den);
    spec.ok = std::isfinite(energy) && spec.residual < 1e-6;
    return spec;
}

inline double interp(const std::vector<double>& x, const std::vector<double>& y, double xq) {
    int n = (int)x.size();
    if (n == 0) return 0.0;
    if (xq <= x.front()) return y.front();
    if (xq >= x.back()) return y.back();
    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (x[mid] <= xq) lo = mid;
        else hi = mid;
    }
    double t = (xq - x[lo]) / (x[hi] - x[lo]);
    return y[lo] * (1.0 - t) + y[hi] * t;
}

}
