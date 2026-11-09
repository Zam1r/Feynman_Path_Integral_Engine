#pragma once
#include "action.hpp"
#include "rng.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace feynman {

struct MetroConfig {
    int N;
    double T;
    double mass;
    double hbar;
    Potential pot;
    int sweeps;
    int therm;
    int tune;
    double step;
    int chains;
    uint64_t seed;
    int bins;
    double xmin;
    double xmax;
    int save_paths;
    int batches;
    double x0;
    bool do_corr;
    bool do_world;
    bool do_trace;
};

struct MetroResult {
    std::vector<double> x;
    std::vector<double> psi2;
    std::vector<double> err;
    std::vector<double> sweep;
    std::vector<double> trace_x2;
    std::vector<double> trace_acc;
    std::vector<double> trace_action;
    std::vector<int> trace_measuring;
    std::vector<double> corr_tau;
    std::vector<double> corr_xx;
    std::vector<double> corr_err;
    std::vector<std::vector<double>> paths;
    std::vector<double> world;
    int world_bins;
    double world_xmin;
    double world_xmax;
    double acceptance;
    double mean_x2;
    double mean_x2_err;
    double mean_energy;
    double mean_energy_err;
    double mean_crossings;
    double step;
    double psi2_at_1;
    double err_at_1;
    double eps;
    int N;
    double T;
    bool ok;
};

inline int bin_of(double x, double xmin, double dx, int bins) {
    return (int)std::floor((x - xmin) / dx);
}

inline double density_at(const std::vector<double>& hist, double samples, double xmin, double dx, int bins, double xq) {
    if (!(samples > 0.0)) return 0.0;
    double pos = (xq - xmin) / dx - 0.5;
    if (pos <= 0.0) return hist[0] / (samples * dx);
    if (pos >= (double)(bins - 1)) return hist[bins - 1] / (samples * dx);
    int i = (int)std::floor(pos);
    double t = pos - (double)i;
    double y0 = hist[i] / (samples * dx);
    double y1 = hist[i + 1] / (samples * dx);
    return y0 * (1.0 - t) + y1 * t;
}

inline MetroResult metropolis_ground_state(const MetroConfig& cfg) {
    MetroResult R;
    R.ok = false;
    R.acceptance = 0.0;
    R.mean_x2 = 0.0;
    R.mean_x2_err = 0.0;
    R.mean_energy = 0.0;
    R.mean_energy_err = 0.0;
    R.mean_crossings = 0.0;
    R.step = cfg.step;
    R.psi2_at_1 = 0.0;
    R.err_at_1 = 0.0;
    R.N = cfg.N;
    R.T = cfg.T;
    R.eps = cfg.T / (double)cfg.N;
    R.world_bins = cfg.bins;
    R.world_xmin = cfg.xmin;
    R.world_xmax = cfg.xmax;
    if (cfg.N < 2 || cfg.sweeps < cfg.batches || cfg.batches < 2 || cfg.bins < 4 || !(cfg.xmax > cfg.xmin)) return R;
    double eps = R.eps;
    double dx = (cfg.xmax - cfg.xmin) / (double)cfg.bins;
    int chains = std::max(1, cfg.chains);
    int nbatch = cfg.batches;
    int per = cfg.sweeps / nbatch;
    if (per < 1) return R;
    std::vector<std::vector<double>> batchHist(nbatch, std::vector<double>(cfg.bins, 0.0));
    std::vector<std::vector<double>> batchCorr(nbatch, std::vector<double>(cfg.N, 0.0));
    std::vector<double> batchX2(nbatch, 0.0);
    std::vector<double> batchE(nbatch, 0.0);
    std::vector<double> batchCross(nbatch, 0.0);
    std::vector<double> batchCount(nbatch, 0.0);
    std::vector<double> corrCount(nbatch, 0.0);
    std::vector<double> world;
    if (cfg.do_world) world.assign((size_t)cfg.N * (size_t)cfg.bins, 0.0);
    double world_count = 0.0;
    double accept_sum = 0.0;
    double trial_sum = 0.0;
    double step = cfg.step;
    if (!(step > 0.0)) step = 2.0 * std::sqrt(eps * cfg.hbar / cfg.mass);
    int trace_stride = 1;
    int trace_budget = cfg.tune + cfg.therm + cfg.sweeps;
    if (trace_budget > 2000) trace_stride = (trace_budget + 1999) / 2000;

    for (int c = 0; c < chains; ++c) {
        Rng rng(cfg.seed + 0x100000001b3ULL * (uint64_t)(c + 1));
        std::vector<double> path(cfg.N, cfg.x0);
        if (cfg.pot.kind == DOUBLEWELL) {
            double side = (c % 2 == 0) ? cfg.pot.b : -cfg.pot.b;
            std::fill(path.begin(), path.end(), side);
        }
        double chain_step = step;
        auto sweep_once = [&](bool record) -> double {
            int acc = 0;
            for (int k = 0; k < cfg.N; ++k) {
                int prev = k - 1;
                if (prev < 0) prev = cfg.N - 1;
                int next = k + 1;
                if (next == cfg.N) next = 0;
                double prop = path[k] + (rng.uniform() * 2.0 - 1.0) * chain_step;
                double oldS = link_action(path[prev], path[k], eps, cfg.mass, cfg.hbar, cfg.pot)
                            + link_action(path[k], path[next], eps, cfg.mass, cfg.hbar, cfg.pot);
                double newS = link_action(path[prev], prop, eps, cfg.mass, cfg.hbar, cfg.pot)
                            + link_action(prop, path[next], eps, cfg.mass, cfg.hbar, cfg.pot);
                double dS = newS - oldS;
                if (dS <= 0.0 || (dS < 40.0 && rng.uniform() < std::exp(-dS))) {
                    path[k] = prop;
                    ++acc;
                }
            }
            if (record) {
                accept_sum += (double)acc;
                trial_sum += (double)cfg.N;
            }
            return (double)acc / (double)cfg.N;
        };
        int tune_acc = 0;
        int tune_try = 0;
        for (int s = 0; s < cfg.tune; ++s) {
            double rate_one = sweep_once(false);
            tune_acc += (int)std::lround(rate_one * cfg.N);
            tune_try += cfg.N;
            if ((s + 1) % 20 == 0 && tune_try > 0) {
                double rate = (double)tune_acc / (double)tune_try;
                if (rate > 0.55) chain_step *= 1.1;
                else if (rate < 0.32) chain_step *= 0.9;
                if (chain_step < 1e-4) chain_step = 1e-4;
                if (chain_step > 6.0) chain_step = 6.0;
                tune_acc = 0;
                tune_try = 0;
            }
            if (cfg.do_trace && c == 0 && (s % trace_stride) == 0) {
                double x2 = 0.0;
                for (double v : path) x2 += v * v;
                R.sweep.push_back((double)s);
                R.trace_x2.push_back(x2 / (double)cfg.N);
                R.trace_acc.push_back(rate_one);
                R.trace_action.push_back(periodic_action(path, eps, cfg.mass, cfg.hbar, cfg.pot));
                R.trace_measuring.push_back(0);
            }
        }
        for (int s = 0; s < cfg.therm; ++s) {
            double rate_one = sweep_once(false);
            if (cfg.do_trace && c == 0 && (s % trace_stride) == 0) {
                double x2 = 0.0;
                for (double v : path) x2 += v * v;
                R.sweep.push_back((double)(cfg.tune + s));
                R.trace_x2.push_back(x2 / (double)cfg.N);
                R.trace_acc.push_back(rate_one);
                R.trace_action.push_back(periodic_action(path, eps, cfg.mass, cfg.hbar, cfg.pot));
                R.trace_measuring.push_back(0);
            }
        }
        int saved = 0;
        int want = std::max(0, cfg.save_paths);
        int save_every = std::max(1, (per * nbatch) / std::max(1, want));
        for (int s = 0; s < per * nbatch; ++s) {
            double rate_one = sweep_once(true);
            int b = s / per;
            if (b >= nbatch) b = nbatch - 1;
            double x2 = 0.0;
            for (int j = 0; j < cfg.N; ++j) {
                double v = path[j];
                x2 += v * v;
                int bin = bin_of(v, cfg.xmin, dx, cfg.bins);
                if (bin >= 0 && bin < cfg.bins) {
                    batchHist[b][bin] += 1.0;
                    if (cfg.do_world && c == 0) world[(size_t)j * cfg.bins + bin] += 1.0;
                }
            }
            x2 /= (double)cfg.N;
            batchX2[b] += x2;
            batchE[b] += energy_estimator(path, eps, cfg.mass, cfg.hbar, cfg.pot);
            batchCross[b] += (double)sign_crossings(path);
            batchCount[b] += 1.0;
            if (cfg.do_corr && c == 0) {
                corrCount[b] += 1.0;
                for (int dt = 0; dt < cfg.N; ++dt) {
                    double acc = 0.0;
                    int k = dt;
                    for (int j = 0; j < cfg.N; ++j) {
                        acc += path[j] * path[k];
                        ++k;
                        if (k == cfg.N) k = 0;
                    }
                    batchCorr[b][dt] += acc / (double)cfg.N;
                }
            }
            if (c == 0) {
                if (cfg.do_world) world_count += 1.0;
                if (cfg.do_trace && (s % trace_stride) == 0) {
                    R.sweep.push_back((double)(cfg.tune + cfg.therm + s));
                    R.trace_x2.push_back(x2);
                    R.trace_acc.push_back(rate_one);
                    R.trace_action.push_back(periodic_action(path, eps, cfg.mass, cfg.hbar, cfg.pot));
                    R.trace_measuring.push_back(1);
                }
                if (want > 0 && saved < want && (s % save_every) == 0) {
                    R.paths.push_back(path);
                    ++saved;
                }
            }
        }
        step = chain_step;
    }

    R.step = step;
    R.acceptance = (trial_sum > 0.0) ? accept_sum / trial_sum : 0.0;
    R.x.resize(cfg.bins);
    R.psi2.assign(cfg.bins, 0.0);
    R.err.assign(cfg.bins, 0.0);
    std::vector<double> batchAt1(nbatch, 0.0);
    std::vector<double> batchMeanX2(nbatch, 0.0);
    std::vector<double> batchMeanE(nbatch, 0.0);
    std::vector<double> batchMeanC(nbatch, 0.0);
    for (int b = 0; b < nbatch; ++b) {
        double samples = batchCount[b] * (double)cfg.N;
        if (!(samples > 0.0)) continue;
        for (int i = 0; i < cfg.bins; ++i) R.psi2[i] += batchHist[b][i] / (samples * dx);
        batchAt1[b] = density_at(batchHist[b], samples, cfg.xmin, dx, cfg.bins, 1.0);
        batchMeanX2[b] = batchX2[b] / batchCount[b];
        batchMeanE[b] = batchE[b] / batchCount[b];
        batchMeanC[b] = batchCross[b] / batchCount[b];
    }
    double invB = 1.0 / (double)nbatch;
    for (int i = 0; i < cfg.bins; ++i) {
        R.x[i] = cfg.xmin + ((double)i + 0.5) * dx;
        R.psi2[i] *= invB;
    }
    double scale = 1.0 / std::sqrt((double)nbatch * (double)(nbatch - 1));
    for (int b = 0; b < nbatch; ++b) {
        double samples = batchCount[b] * (double)cfg.N;
        if (!(samples > 0.0)) continue;
        for (int i = 0; i < cfg.bins; ++i) {
            double rho = batchHist[b][i] / (samples * dx);
            double dlt = rho - R.psi2[i];
            R.err[i] += dlt * dlt;
        }
    }
    for (int i = 0; i < cfg.bins; ++i) R.err[i] = std::sqrt(R.err[i]) * scale;
    double mean_at1 = 0.0;
    for (double v : batchAt1) mean_at1 += v;
    mean_at1 *= invB;
    double var_at1 = 0.0;
    for (double v : batchAt1) {
        double dlt = v - mean_at1;
        var_at1 += dlt * dlt;
    }
    R.psi2_at_1 = mean_at1;
    R.err_at_1 = std::sqrt(var_at1) * scale;
    double meanx2 = 0.0;
    for (double v : batchMeanX2) meanx2 += v;
    meanx2 *= invB;
    double varx2 = 0.0;
    for (double v : batchMeanX2) {
        double dlt = v - meanx2;
        varx2 += dlt * dlt;
    }
    R.mean_x2 = meanx2;
    R.mean_x2_err = std::sqrt(varx2) * scale;
    double meane = 0.0;
    for (double v : batchMeanE) meane += v;
    meane *= invB;
    double vare = 0.0;
    for (double v : batchMeanE) {
        double dlt = v - meane;
        vare += dlt * dlt;
    }
    R.mean_energy = meane;
    R.mean_energy_err = std::sqrt(vare) * scale;
    double meanc = 0.0;
    for (double v : batchMeanC) meanc += v;
    R.mean_crossings = meanc * invB;

    if (cfg.do_corr) {
        R.corr_tau.resize(cfg.N);
        R.corr_xx.assign(cfg.N, 0.0);
        R.corr_err.assign(cfg.N, 0.0);
        for (int b = 0; b < nbatch; ++b) {
            double invc = (corrCount[b] > 0.0) ? 1.0 / corrCount[b] : 0.0;
            for (int dt = 0; dt < cfg.N; ++dt) R.corr_xx[dt] += batchCorr[b][dt] * invc;
        }
        for (int dt = 0; dt < cfg.N; ++dt) {
            R.corr_tau[dt] = (double)dt * eps;
            R.corr_xx[dt] *= invB;
        }
        for (int b = 0; b < nbatch; ++b) {
            double invc = (corrCount[b] > 0.0) ? 1.0 / corrCount[b] : 0.0;
            for (int dt = 0; dt < cfg.N; ++dt) {
                double v = batchCorr[b][dt] * invc;
                double dlt = v - R.corr_xx[dt];
                R.corr_err[dt] += dlt * dlt;
            }
        }
        for (int dt = 0; dt < cfg.N; ++dt) R.corr_err[dt] = std::sqrt(R.corr_err[dt]) * scale;
    }
    if (cfg.do_world && world_count > 0.0) {
        double norm = world_count * dx;
        for (double& v : world) v /= norm;
        R.world = std::move(world);
    }
    R.ok = true;
    return R;
}

}
