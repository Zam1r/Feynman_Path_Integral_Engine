#pragma once
#include "action.hpp"
#include "rng.hpp"
#include <vector>
#include <atomic>
#include <thread>
#include <cmath>

namespace feynman {

struct DirectConfig {
    int N;
    double T;
    double mass;
    double hbar;
    Potential pot;
    double xmin;
    double xmax;
    double dx;
    int samples;
    int batches;
    uint64_t seed;
    int threads;
};

struct DirectResult {
    std::vector<double> x;
    std::vector<double> psi2;
    std::vector<double> err;
    double psi2_at_1;
    double err_at_1;
    bool ok;
};

struct BridgePath {
    std::vector<double> x;
    double weight;
};

inline double interp_linear(const std::vector<double>& x, const std::vector<double>& y, double xq) {
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

inline double bridge_potential_weight(double xb, int N, double eps, double mass, double hbar, const Potential& V, Rng& rng, std::vector<double>* path) {
    double sigma2 = eps * hbar / mass;
    double x_prev = xb;
    double spot = 0.0;
    if (path) {
        path->clear();
        path->push_back(xb);
    }
    for (int j = 1; j < N; ++j) {
        double steps = (double)(N - j + 1);
        double mean = x_prev + (xb - x_prev) / steps;
        double var = sigma2 * (steps - 1.0) / steps;
        double x = mean + std::sqrt(var) * rng.normal();
        double mid = 0.5 * (x + x_prev);
        spot += (eps / hbar) * V(mid);
        x_prev = x;
        if (path) path->push_back(x);
    }
    double mid = 0.5 * (xb + x_prev);
    spot += (eps / hbar) * V(mid);
    if (path) path->push_back(xb);
    if (spot > 700.0) return 0.0;
    return std::exp(-spot);
}

inline std::vector<BridgePath> sample_bridges(double xb, int N, double T, double mass, double hbar, const Potential& V, int npaths, uint64_t seed) {
    double eps = T / (double)N;
    Rng rng(seed);
    std::vector<BridgePath> out;
    out.reserve(npaths);
    for (int i = 0; i < npaths; ++i) {
        BridgePath p;
        p.weight = bridge_potential_weight(xb, N, eps, mass, hbar, V, rng, &p.x);
        out.push_back(std::move(p));
    }
    return out;
}

inline DirectResult direct_ground_density(const DirectConfig& cfg) {
    DirectResult result;
    result.ok = false;
    result.psi2_at_1 = 0.0;
    result.err_at_1 = 0.0;
    if (cfg.N < 2 || cfg.samples < cfg.batches || cfg.batches < 2 || !(cfg.dx > 0.0) || !(cfg.xmax > cfg.xmin)) return result;
    int nx = (int)std::floor((cfg.xmax - cfg.xmin) / cfg.dx + 0.5) + 1;
    if (nx < 3) return result;
    std::vector<double> xs(nx);
    for (int i = 0; i < nx; ++i) xs[i] = cfg.xmin + (double)i * cfg.dx;
    int nbatch = cfg.batches;
    int per = cfg.samples / nbatch;
    if (per < 1) return result;
    std::vector<std::vector<double>> batchG(nbatch, std::vector<double>(nx, 0.0));
    double eps = cfg.T / (double)cfg.N;
    int threads = cfg.threads;
    if (threads < 1) threads = 1;
    if (threads > nx) threads = nx;
    std::atomic<int> cursor{0};
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            (void)t;
            while (true) {
                int i = cursor.fetch_add(1);
                if (i >= nx) break;
                Rng rng(cfg.seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(i + 1));
                double xb = xs[i];
                for (int b = 0; b < nbatch; ++b) {
                    double sum = 0.0;
                    for (int s = 0; s < per; ++s) {
                        sum += bridge_potential_weight(xb, cfg.N, eps, cfg.mass, cfg.hbar, cfg.pot, rng, nullptr);
                    }
                    batchG[b][i] = sum / (double)per;
                }
            }
        });
    }
    for (auto& th : pool) th.join();
    result.x = xs;
    result.psi2.assign(nx, 0.0);
    result.err.assign(nx, 0.0);
    std::vector<double> meanG(nx, 0.0);
    for (int b = 0; b < nbatch; ++b) {
        double Z = 0.0;
        Z += 0.5 * batchG[b][0];
        Z += 0.5 * batchG[b][nx - 1];
        for (int i = 1; i < nx - 1; ++i) Z += batchG[b][i];
        Z *= cfg.dx;
        if (!(Z > 0.0)) continue;
        for (int i = 0; i < nx; ++i) {
            double rho = batchG[b][i] / Z;
            result.psi2[i] += rho;
            meanG[i] += batchG[b][i];
        }
    }
    double invB = 1.0 / (double)nbatch;
    for (int i = 0; i < nx; ++i) result.psi2[i] *= invB;
    for (int b = 0; b < nbatch; ++b) {
        double Z = 0.0;
        Z += 0.5 * batchG[b][0];
        Z += 0.5 * batchG[b][nx - 1];
        for (int i = 1; i < nx - 1; ++i) Z += batchG[b][i];
        Z *= cfg.dx;
        if (!(Z > 0.0)) continue;
        for (int i = 0; i < nx; ++i) {
            double rho = batchG[b][i] / Z;
            double d = rho - result.psi2[i];
            result.err[i] += d * d;
        }
    }
    double scale = 1.0 / std::sqrt((double)nbatch * (double)(nbatch - 1));
    for (int i = 0; i < nx; ++i) result.err[i] = std::sqrt(result.err[i]) * scale;
    result.psi2_at_1 = interp_linear(result.x, result.psi2, 1.0);
    result.err_at_1 = interp_linear(result.x, result.err, 1.0);
    result.ok = true;
    return result;
}

}
