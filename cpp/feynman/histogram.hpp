#pragma once
#include "potential.hpp"
#include "rng.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace feynman {

struct LargeHist {
    std::vector<double> x;
    std::vector<double> psi2;
    std::vector<double> err;
    std::vector<double> count;
    double acceptance = 0;
    int N = 0;
    double T = 0;
    int sweeps = 0;
    bool ok = false;
};

inline bool write_histogram_csv(const char* path, const LargeHist& H, const char* pot, double mass, double hbar) {
    if (!H.ok || H.x.empty()) return false;
    FILE* f = std::fopen(path, "w");
    if (!f) return false;
    std::fprintf(f, "N,%d\nT,%.8g\nmass,%.8g\nhbar,%.8g\npotential,%s\nsweeps,%d\nacceptance,%.6f\n", H.N, H.T, mass, hbar, pot ? pot : "ho", H.sweeps, H.acceptance);
    std::fprintf(f, "x,count,density,error\n");
    for (size_t i = 0; i < H.x.size(); ++i) {
        double c = i < H.count.size() ? H.count[i] : 0.0;
        double e = i < H.err.size() ? H.err[i] : 0.0;
        std::fprintf(f, "%.8g,%.8g,%.8g,%.8g\n", H.x[i], c, H.psi2[i], e);
    }
    std::fclose(f);
    return true;
}

inline LargeHist large_histogram(const Potential& V, int N, double T, double mass, double hbar, int sweeps, int therm, int bins, double xmin, double xmax, uint64_t seed) {
    LargeHist H;
    H.acceptance = 0.0;
    H.N = N;
    H.T = T;
    H.sweeps = sweeps;
    H.ok = false;
    if (N < 8 || sweeps < 4 || bins < 8 || !(T > 0.0) || !(mass > 0.0) || !(hbar > 0.0) || !(xmax > xmin)) return H;
    (void)therm;
    double eps = T / (double)N;
    double dx = (xmax - xmin) / (double)bins;
    double cv = eps / hbar;
    const double pi = 3.14159265358979323846;
    int kmax = 48;
    if (kmax > N / 2) kmax = N / 2;
    int nbatch = 4;
    if (sweeps < nbatch) nbatch = 2;
    int per = sweeps / nbatch;
    if (per < 1) return H;
    sweeps = per * nbatch;
    H.sweeps = sweeps;
    int shifts = 24;
    Rng rng(seed);
    std::vector<double> path;
    std::vector<std::vector<double>> batch(nbatch, std::vector<double>(bins, 0.0));
    double shift_w = 0.65;
    int acc = 0;
    int trials = 0;
    int stride = 1;
    if (N > 4096) stride = N / 4096;
    int nkeep = N / stride;
    if (nkeep < 8) return H;
    path.assign(nkeep, 0.0);
    auto draw_free = [&]() {
        std::fill(path.begin(), path.end(), 0.0);
        for (int k = 1; k <= kmax; ++k) {
            double s = std::sin(pi * (double)k / (double)N);
            double alpha_kin = (mass * (double)N * (double)N * s * s) / (T * hbar);
            double alpha_pot = 0.0;
            if (V.kind == HO) alpha_pot = (T * V.a) / (4.0 * hbar);
            double alpha = alpha_kin + alpha_pot;
            if (!(alpha > 0.0)) continue;
            double var = 1.0 / (2.0 * alpha);
            double sig = std::sqrt(var);
            double A = sig * rng.normal();
            double B = sig * rng.normal();
            double w = 2.0 * pi * (double)k / (double)N;
            for (int i = 0; i < nkeep; ++i) {
                double th = w * (double)(i * stride);
                path[i] += A * std::cos(th) + B * std::sin(th);
            }
        }
        if (V.kind == HO) {
            double sig0 = std::sqrt(hbar / (V.a * T));
            double a0 = sig0 * rng.normal();
            for (int i = 0; i < nkeep; ++i) path[i] += a0;
        } else if (V.kind == DOUBLEWELL) {
            double side = ((rng.next_u64() & 1ULL) == 0ULL) ? V.b : -V.b;
            for (int i = 0; i < nkeep; ++i) path[i] += side;
        }
    };
    auto shift_once = [&](bool record) {
        double delta = (rng.uniform() * 2.0 - 1.0) * shift_w;
        double dS = 0.0;
        double weight = cv * (double)stride;
        for (int i = 0; i < nkeep; ++i) dS += weight * (V(path[i] + delta) - V(path[i]));
        bool take = dS <= 0.0 || (dS < 80.0 && rng.uniform() < std::exp(-dS));
        if (take) {
            for (int i = 0; i < nkeep; ++i) path[i] += delta;
        }
        if (record) {
            trials++;
            if (take) acc++;
        }
        return take;
    };
    for (int s = 0; s < sweeps; ++s) {
        draw_free();
        int got = 0;
        if (V.kind != HO) {
            for (int t = 0; t < shifts; ++t) {
                if (shift_once(true)) got++;
            }
        }
        if (s < 4) {
            double rate = (double)got / (double)shifts;
            if (rate > 0.55) shift_w *= 1.1;
            else if (rate < 0.25) shift_w *= 0.85;
            if (shift_w < 0.05) shift_w = 0.05;
            if (shift_w > 2.5) shift_w = 2.5;
        }
        int b = s / per;
        if (b >= nbatch) b = nbatch - 1;
        for (int i = 0; i < nkeep; ++i) {
            int bin = (int)((path[i] - xmin) / dx);
            if (bin >= 0 && bin < bins) batch[b][bin] += 1.0;
        }
    }
    H.acceptance = (trials > 0) ? (double)acc / (double)trials : 1.0;
    H.x.resize(bins);
    H.psi2.assign(bins, 0.0);
    H.err.assign(bins, 0.0);
    H.count.assign(bins, 0.0);
    std::vector<double> mean(bins, 0.0);
    double samples = (double)per * (double)nkeep;
    for (int b = 0; b < nbatch; ++b) {
        for (int i = 0; i < bins; ++i) {
            H.count[i] += batch[b][i];
            mean[i] += batch[b][i] / (samples * dx);
        }
    }
    double invB = 1.0 / (double)nbatch;
    for (int i = 0; i < bins; ++i) {
        H.x[i] = xmin + ((double)i + 0.5) * dx;
        H.psi2[i] = mean[i] * invB;
    }
    double scale = 1.0 / std::sqrt((double)nbatch * (double)(nbatch - 1));
    for (int b = 0; b < nbatch; ++b) {
        for (int i = 0; i < bins; ++i) {
            double rho = batch[b][i] / (samples * dx);
            double d = rho - H.psi2[i];
            H.err[i] += d * d;
        }
    }
    for (int i = 0; i < bins; ++i) H.err[i] = std::sqrt(H.err[i]) * scale;
    H.ok = true;
    return H;
}

}
