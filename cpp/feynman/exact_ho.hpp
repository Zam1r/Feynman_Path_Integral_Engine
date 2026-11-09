#pragma once
#include <vector>
#include <cmath>

namespace feynman {

struct HoLattice {
    double alpha;
    bool ok;
};

inline HoLattice lattice_alpha(int N, double T, double mass, double hbar, double k) {
    HoLattice out{0.0, false};
    if (N < 1 || T <= 0.0 || mass <= 0.0 || hbar <= 0.0 || k <= 0.0) return out;
    double eps = T / (double)N;
    double akin = mass / (2.0 * eps * hbar);
    double bpot = (eps * k) / (8.0 * hbar);
    double diag = akin + bpot;
    double cross = -2.0 * akin + 2.0 * bpot;
    int n = N - 1;
    if (n == 0) {
        out.alpha = 4.0 * bpot;
        out.ok = out.alpha > 0.0;
        return out;
    }
    std::vector<double> M((size_t)n * (size_t)n, 0.0);
    std::vector<double> J(n, 0.0);
    double cnst = 0.0;
    for (int link = 1; link <= N; ++link) {
        int left = link - 1;
        int right = link;
        bool lb = left == 0;
        bool rb = right == N;
        if (!lb) M[(size_t)(left - 1) * n + (left - 1)] += diag;
        else cnst += diag;
        if (!rb) M[(size_t)(right - 1) * n + (right - 1)] += diag;
        else cnst += diag;
        if (!lb && !rb) {
            int a = left - 1;
            int b = right - 1;
            M[(size_t)a * n + b] += 0.5 * cross;
            M[(size_t)b * n + a] += 0.5 * cross;
        } else if (lb && !rb) {
            J[right - 1] += cross;
        } else if (!lb && rb) {
            J[left - 1] += cross;
        } else {
            cnst += cross;
        }
    }
    std::vector<double> L = M;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = L[(size_t)i * n + j];
            for (int p = 0; p < j; ++p) sum -= L[(size_t)i * n + p] * L[(size_t)j * n + p];
            if (i == j) {
                if (sum <= 1e-18) return out;
                L[(size_t)i * n + j] = std::sqrt(sum);
            } else {
                L[(size_t)i * n + j] = sum / L[(size_t)j * n + j];
            }
        }
    }
    std::vector<double> w(n), z(n);
    for (int i = 0; i < n; ++i) {
        double sum = J[i];
        for (int p = 0; p < i; ++p) sum -= L[(size_t)i * n + p] * w[p];
        w[i] = sum / L[(size_t)i * n + i];
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = w[i];
        for (int p = i + 1; p < n; ++p) sum -= L[(size_t)p * n + i] * z[p];
        z[i] = sum / L[(size_t)i * n + i];
    }
    double quad = 0.0;
    for (int i = 0; i < n; ++i) quad += J[i] * z[i];
    out.alpha = cnst - 0.25 * quad;
    out.ok = out.alpha > 0.0 && std::isfinite(out.alpha);
    return out;
}

inline double gaussian_psi2(double x, double alpha) {
    return std::sqrt(alpha / 3.14159265358979323846) * std::exp(-alpha * x * x);
}

}
