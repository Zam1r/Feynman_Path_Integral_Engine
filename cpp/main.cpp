#include "feynman/action.hpp"
#include "feynman/direct.hpp"
#include "feynman/exact_ho.hpp"
#include "feynman/metropolis.hpp"
#include "feynman/potential.hpp"
#include "feynman/rng.hpp"
#include "feynman/schrodinger.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace feynman;

static void usage() {
    std::cerr << "usage:\n"
              << "  feynman selftest\n"
              << "  feynman suite <outdir>\n"
              << "  feynman run --pot ho --N 36 --T 8 --sweeps 8000 --mass 1 --hbar 1\n"
              << "             --chains 2 --seed 1 --xmin -6 --xmax 6 --bins 140\n";
}

static std::string arg(int argc, char** argv, const std::string& key, const std::string& fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (key == argv[i]) return argv[i + 1];
    }
    return fallback;
}

static double argd(int argc, char** argv, const std::string& key, double fallback) {
    std::string s = arg(argc, argv, key, "");
    if (s.empty()) return fallback;
    return std::atof(s.c_str());
}

static int argi(int argc, char** argv, const std::string& key, int fallback) {
    std::string s = arg(argc, argv, key, "");
    if (s.empty()) return fallback;
    return std::atoi(s.c_str());
}

static MetroConfig base_cfg() {
    MetroConfig c{};
    c.N = 32;
    c.T = 8.0;
    c.mass = 1.0;
    c.hbar = 1.0;
    c.pot = parse_potential("ho");
    c.sweeps = 8000;
    c.therm = 1000;
    c.tune = 400;
    c.step = 0.0;
    c.chains = 2;
    c.seed = 1;
    c.bins = 160;
    c.xmin = -6.0;
    c.xmax = 6.0;
    c.save_paths = 0;
    c.batches = 8;
    c.x0 = 0.0;
    c.do_corr = false;
    c.do_world = false;
    c.do_trace = false;
    return c;
}

static bool open_out(std::ofstream& out, const fs::path& path) {
    out.open(path);
    if (!out) {
        std::cerr << "cannot write " << path.string() << "\n";
        return false;
    }
    out.setf(std::ios::fmtflags(0), std::ios::floatfield);
    out.precision(12);
    return true;
}

static int selftest() {
    int fails = 0;
    auto check = [&](const std::string& name, bool ok, double got, double want) {
        double err = std::abs(got - want);
        std::cout << (ok ? "ok  " : "FAIL") << "  " << name << "  got " << got << "  want " << want << "\n";
        if (!ok) ++fails;
        (void)err;
    };

    std::vector<double> d{2.0, 2.0, 2.0};
    std::vector<double> off{-1.0, -1.0};
    double e0 = 2.0 - std::sqrt(2.0);
    check("sturm below", sturm_count(d, off, e0 - 1e-8) == 0, (double)sturm_count(d, off, e0 - 1e-8), 0.0);
    check("sturm above", sturm_count(d, off, e0 + 1e-6) == 1, (double)sturm_count(d, off, e0 + 1e-6), 1.0);

    Potential ho = parse_potential("ho");
    Spectrum sp = solve_schrodinger(ho, 1.0, 1.0, -12.0, 12.0, 1600);
    double analytic1 = ho_ground_psi2(1.0, 1.0, 1.0, 1.0);
    double se1 = sp.ok ? interp(sp.x, sp.psi2, 1.0) : 0.0;
    std::cout << "schrodinger residual " << sp.residual << " energy " << sp.energy << "\n";
    check("ho energy", sp.ok && std::abs(sp.energy - 0.5) < 2e-3, sp.energy, 0.5);
    check("ho psi2(1)", sp.ok && std::abs(se1 - analytic1) < 2e-3, se1, analytic1);

    HoLattice fine = lattice_alpha(160, 24.0, 1.0, 1.0, 1.0);
    check("lattice alpha large N", fine.ok && std::abs(fine.alpha - 1.0) < 0.02, fine.alpha, 1.0);
    HoLattice mile = lattice_alpha(7, 4.0, 1.0, 1.0, 1.0);
    double mile_psi = mile.ok ? gaussian_psi2(1.0, mile.alpha) : 0.0;
    std::cout << "milestone lattice alpha " << mile.alpha << " psi2(1) " << mile_psi << "\n";
    check("milestone lattice finite", mile.ok && mile_psi > 0.12 && mile_psi < 0.35, mile_psi, analytic1);

    auto bridges = sample_bridges(0.0, 80, 4.0, 1.0, 1.0, ho, 3000, 7);
    double mean = 0.0;
    for (const auto& p : bridges) mean += p.x[40];
    mean /= (double)bridges.size();
    double var = 0.0;
    for (const auto& p : bridges) {
        double dlt = p.x[40] - mean;
        var += dlt * dlt;
    }
    var /= (double)bridges.size();
    check("bridge variance", std::abs(var - 1.0) < 0.08 && std::abs(mean) < 0.05, var, 1.0);

    DirectConfig dc{};
    dc.N = 7;
    dc.T = 4.0;
    dc.mass = 1.0;
    dc.hbar = 1.0;
    dc.pot = ho;
    dc.xmin = -2.0;
    dc.xmax = 2.0;
    dc.dx = 0.5;
    dc.samples = 8000;
    dc.batches = 8;
    dc.seed = 11;
    dc.threads = 2;
    DirectResult dr = direct_ground_density(dc);
    double exact_at = gaussian_psi2(0.0, mile.alpha);
    check("direct matches lattice at 0", dr.ok && std::abs(dr.psi2[4] - exact_at) < 0.03, dr.ok ? dr.psi2[4] : 0.0, exact_at);

    if (fails) std::cout << fails << " checks failed\n";
    else std::cout << "all checks passed\n";
    return fails ? 1 : 0;
}

static double fit_gap(const MetroResult& R) {
    double n = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (size_t i = 0; i < R.corr_tau.size(); ++i) {
        double tau = R.corr_tau[i];
        double c = R.corr_xx[i];
        if (tau < 0.35 || tau > 2.4 || !(c > 0.0)) continue;
        double y = std::log(c);
        sx += tau;
        sy += y;
        sxx += tau * tau;
        sxy += tau * y;
        n += 1.0;
    }
    if (n < 3.0) return 0.0;
    double den = n * sxx - sx * sx;
    if (!(std::abs(den) > 0.0)) return 0.0;
    double slope = (n * sxy - sx * sy) / den;
    return -slope;
}

static double tail_prob(const std::vector<double>& x, const std::vector<double>& y, double edge) {
    if (x.size() < 2) return 0.0;
    double p = 0.0;
    for (size_t i = 0; i + 1 < x.size(); ++i) {
        double mid = 0.5 * (x[i] + x[i + 1]);
        double w = x[i + 1] - x[i];
        if (std::abs(mid) > edge) p += 0.5 * (y[i] + y[i + 1]) * w;
    }
    return p;
}

static void write_density_row(std::ofstream& out, double x, double psi2, double err, double exact, double analytic, double finiteT) {
    out << x << "," << psi2 << "," << err << "," << exact << "," << analytic << "," << finiteT << "\n";
}

static int suite(const fs::path& dir) {
    auto t0 = std::chrono::steady_clock::now();
    if (selftest() != 0) return 1;
    fs::create_directories(dir);
    Potential ho = parse_potential("ho");
    const double analytic1 = ho_ground_psi2(1.0, 1.0, 1.0, 1.0);
    const double finite1 = ho_finite_T_psi2(1.0, 1.0, 1.0, 1.0, 4.0);

    std::cout << "schrodinger...\n";
    Spectrum sp = solve_schrodinger(ho, 1.0, 1.0, -10.0, 10.0, 4001);
    if (!sp.ok) {
        std::cerr << "schrodinger failed, residual " << sp.residual << "\n";
        return 1;
    }
    {
        std::ofstream out;
        if (!open_out(out, dir / "schrodinger_ho.csv")) return 1;
        out << "x,psi,psi2,analytic\n";
        for (size_t i = 0; i < sp.x.size(); ++i) {
            out << sp.x[i] << "," << sp.psi[i] << "," << sp.psi2[i] << "," << ho_ground_psi2(sp.x[i], 1.0, 1.0, 1.0) << "\n";
        }
    }
    double se1 = interp(sp.x, sp.psi2, 1.0);

    std::cout << "exact lattice family...\n";
    HoLattice mile = lattice_alpha(7, 4.0, 1.0, 1.0, 1.0);
    {
        std::ofstream out;
        if (!open_out(out, dir / "exact_alpha.csv")) return 1;
        out << "N,T,alpha,psi2_1,ground_1,finiteT_1\n";
        int Ns[] = {4, 7, 12, 24, 48, 96, 160};
        double Ts[] = {2.0, 3.0, 4.0, 8.0, 12.0, 16.0};
        for (int N : Ns) {
            for (double T : Ts) {
                HoLattice a = lattice_alpha(N, T, 1.0, 1.0, 1.0);
                if (!a.ok) continue;
                out << N << "," << T << "," << a.alpha << "," << gaussian_psi2(1.0, a.alpha) << ","
                    << analytic1 << "," << ho_finite_T_psi2(1.0, 1.0, 1.0, 1.0, T) << "\n";
            }
        }
    }

    std::cout << "direct monte carlo N=7 T=4...\n";
    DirectConfig dc{};
    dc.N = 7;
    dc.T = 4.0;
    dc.mass = 1.0;
    dc.hbar = 1.0;
    dc.pot = ho;
    dc.xmin = -5.5;
    dc.xmax = 5.5;
    dc.dx = 0.05;
    dc.samples = 160000;
    dc.batches = 16;
    dc.seed = 42;
    dc.threads = 4;
    DirectResult dr = direct_ground_density(dc);
    if (!dr.ok) {
        std::cerr << "direct failed\n";
        return 1;
    }
    {
        std::ofstream out;
        if (!open_out(out, dir / "direct_milestone.csv")) return 1;
        out << "x,psi2,err,exact,analytic,finiteT\n";
        for (size_t i = 0; i < dr.x.size(); ++i) {
            write_density_row(out, dr.x[i], dr.psi2[i], dr.err[i], gaussian_psi2(dr.x[i], mile.alpha),
                              ho_ground_psi2(dr.x[i], 1.0, 1.0, 1.0), ho_finite_T_psi2(dr.x[i], 1.0, 1.0, 1.0, 4.0));
        }
    }
    {
        auto bridges = sample_bridges(1.0, 7, 4.0, 1.0, 1.0, ho, 36, 99);
        std::ofstream out;
        if (!open_out(out, dir / "bridges.csv")) return 1;
        out << "path,tau,x,weight\n";
        double eps = 4.0 / 7.0;
        for (size_t p = 0; p < bridges.size(); ++p) {
            for (size_t j = 0; j < bridges[p].x.size(); ++j) {
                out << p << "," << (double)j * eps << "," << bridges[p].x[j] << "," << bridges[p].weight << "\n";
            }
        }
    }

    std::cout << "metropolis N=7...\n";
    MetroConfig small = base_cfg();
    small.N = 7;
    small.T = 4.0;
    small.sweeps = 24000;
    small.therm = 2000;
    small.tune = 600;
    small.chains = 4;
    small.batches = 12;
    small.bins = 220;
    small.xmin = -5.5;
    small.xmax = 5.5;
    small.seed = 3;
    MetroResult ms = metropolis_ground_state(small);
    {
        std::ofstream out;
        if (!open_out(out, dir / "metro_n7.csv")) return 1;
        out << "x,psi2,err\n";
        for (size_t i = 0; i < ms.x.size(); ++i) out << ms.x[i] << "," << ms.psi2[i] << "," << ms.err[i] << "\n";
    }

    std::cout << "metropolis production...\n";
    MetroConfig prod = base_cfg();
    prod.N = 64;
    prod.T = 12.0;
    prod.sweeps = 16000;
    prod.therm = 2000;
    prod.tune = 600;
    prod.chains = 4;
    prod.batches = 10;
    prod.bins = 240;
    prod.xmin = -6.0;
    prod.xmax = 6.0;
    prod.seed = 5;
    prod.save_paths = 16;
    prod.do_corr = true;
    prod.do_world = true;
    prod.do_trace = true;
    MetroResult mp = metropolis_ground_state(prod);
    {
        std::ofstream out;
        if (!open_out(out, dir / "metro_ho_hist.csv")) return 1;
        out << "x,psi2,err\n";
        for (size_t i = 0; i < mp.x.size(); ++i) out << mp.x[i] << "," << mp.psi2[i] << "," << mp.err[i] << "\n";
    }
    {
        std::ofstream out;
        if (!open_out(out, dir / "metro_ho_trace.csv")) return 1;
        out << "sweep,x2,acceptance,action,measuring\n";
        for (size_t i = 0; i < mp.sweep.size(); ++i) {
            out << mp.sweep[i] << "," << mp.trace_x2[i] << "," << mp.trace_acc[i] << "," << mp.trace_action[i] << ","
                << mp.trace_measuring[i] << "\n";
        }
    }
    {
        std::ofstream out;
        if (!open_out(out, dir / "metro_ho_paths.csv")) return 1;
        out << "path,tau,x\n";
        for (size_t p = 0; p < mp.paths.size(); ++p) {
            for (size_t j = 0; j < mp.paths[p].size(); ++j) {
                out << p << "," << (double)j * mp.eps << "," << mp.paths[p][j] << "\n";
            }
        }
    }
    {
        std::ofstream out;
        if (!open_out(out, dir / "metro_ho_corr.csv")) return 1;
        out << "tau,xx,err,theory\n";
        for (size_t i = 0; i < mp.corr_tau.size(); ++i) {
            double tau = mp.corr_tau[i];
            double theory = 0.5 * std::cosh(0.5 * mp.T - tau) / std::sinh(0.5 * mp.T);
            out << tau << "," << mp.corr_xx[i] << "," << mp.corr_err[i] << "," << theory << "\n";
        }
    }
    {
        std::ofstream out;
        if (!open_out(out, dir / "metro_ho_world.csv")) return 1;
        out << "tau,x,density\n";
        double dx = (mp.world_xmax - mp.world_xmin) / (double)mp.world_bins;
        for (int j = 0; j < mp.N; ++j) {
            for (int b = 0; b < mp.world_bins; ++b) {
                double dens = mp.world[(size_t)j * mp.world_bins + b];
                if (dens <= 0.0) continue;
                double x = mp.world_xmin + ((double)b + 0.5) * dx;
                out << (double)j * mp.eps << "," << x << "," << dens << "\n";
            }
        }
    }

    std::cout << "hbar scan...\n";
    {
        std::ofstream out;
        if (!open_out(out, dir / "hbar_scan.csv")) return 1;
        out << "hbar,x,psi2,err,analytic\n";
        double hbars[] = {0.4, 0.7, 1.0, 1.6, 2.4, 3.5};
        for (double hb : hbars) {
            MetroConfig c = base_cfg();
            c.N = 40;
            c.T = 10.0;
            c.hbar = hb;
            c.sweeps = 7000;
            c.therm = 800;
            c.tune = 400;
            c.chains = 2;
            c.batches = 8;
            c.bins = 180;
            c.xmin = -7.0;
            c.xmax = 7.0;
            c.seed = 100 + (uint64_t)std::lround(hb * 10.0);
            std::cout << "  hbar " << hb << "\n";
            MetroResult r = metropolis_ground_state(c);
            for (size_t i = 0; i < r.x.size(); ++i) {
                out << hb << "," << r.x[i] << "," << r.psi2[i] << "," << r.err[i] << ","
                    << ho_ground_psi2(r.x[i], 1.0, hb, 1.0) << "\n";
            }
        }
    }

    std::cout << "other potentials...\n";
    {
        std::ofstream out;
        if (!open_out(out, dir / "potentials.csv")) return 1;
        out << "name,x,V,psi2,err,se_psi2,energy\n";
        struct Job {
            const char* name;
            double xmin, xmax;
        };
        Job jobs[] = {{"anharm", -3.6, 3.6}, {"doublewell", -3.4, 3.4}, {"morse", -1.6, 6.2}};
        uint64_t seed = 200;
        for (const Job& job : jobs) {
            Potential pot = parse_potential(job.name);
            std::cout << "  " << job.name << "\n";
            Spectrum loc = solve_schrodinger(pot, 1.0, 1.0, job.xmin, job.xmax, 1801);
            MetroConfig c = base_cfg();
            c.pot = pot;
            c.N = 48;
            c.T = 10.0;
            c.sweeps = 12000;
            c.therm = 1500;
            c.tune = 700;
            c.chains = 4;
            c.batches = 8;
            c.bins = 180;
            c.xmin = job.xmin;
            c.xmax = job.xmax;
            c.seed = seed;
            seed += 17;
            MetroResult r = metropolis_ground_state(c);
            for (size_t i = 0; i < r.x.size(); ++i) {
                double se = loc.ok ? interp(loc.x, loc.psi2, r.x[i]) : 0.0;
                out << job.name << "," << r.x[i] << "," << pot(r.x[i]) << "," << r.psi2[i] << "," << r.err[i] << ","
                    << se << "," << loc.energy << "\n";
            }
        }
    }

    double gap = fit_gap(mp);
    double forbid = tail_prob(mp.x, mp.psi2, 1.0);
    double dx_small = 0.0;
    if (ms.x.size() > 2) dx_small = ms.x[1] - ms.x[0];
    (void)dx_small;

    {
        std::ofstream out;
        if (!open_out(out, dir / "summary.csv")) return 1;
        out << "key,value\n";
        auto kv = [&](const std::string& k, double v) { out << k << "," << v << "\n"; };
        kv("analytic_psi2_1", analytic1);
        kv("finiteT_psi2_1", finite1);
        kv("schrodinger_energy", sp.energy);
        kv("schrodinger_residual", sp.residual);
        kv("schrodinger_psi2_1", se1);
        kv("lattice_alpha_N7_T4", mile.alpha);
        kv("lattice_psi2_1", gaussian_psi2(1.0, mile.alpha));
        kv("direct_psi2_1", dr.psi2_at_1);
        kv("direct_psi2_1_err", dr.err_at_1);
        kv("direct_N", 7);
        kv("direct_T", 4);
        kv("metro_n7_psi2_1", ms.psi2_at_1);
        kv("metro_n7_psi2_1_err", ms.err_at_1);
        kv("metro_n7_x2", ms.mean_x2);
        kv("metro_n7_accept", ms.acceptance);
        kv("metro_psi2_1", mp.psi2_at_1);
        kv("metro_psi2_1_err", mp.err_at_1);
        kv("metro_x2", mp.mean_x2);
        kv("metro_x2_err", mp.mean_x2_err);
        kv("metro_accept", mp.acceptance);
        kv("metro_step", mp.step);
        kv("metro_N", (double)mp.N);
        kv("metro_T", mp.T);
        kv("metro_eps", mp.eps);
        kv("corr_gap", gap);
        kv("forbidden_prob", forbid);
        kv("forbidden_analytic", std::erfc(1.0));
        kv("ground_x2", 0.5);
    }

    auto t1 = std::chrono::steady_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "analytic |phi(1)|^2 = " << analytic1 << "\n";
    std::cout << "schrodinger |phi(1)|^2 = " << se1 << "   E0 = " << sp.energy << "\n";
    std::cout << "lattice N=7 T=4 |phi(1)|^2 = " << gaussian_psi2(1.0, mile.alpha) << "\n";
    std::cout << "direct |phi(1)|^2 = " << dr.psi2_at_1 << " +/- " << dr.err_at_1 << "\n";
    std::cout << "metro N=7 |phi(1)|^2 = " << ms.psi2_at_1 << " +/- " << ms.err_at_1 << "\n";
    std::cout << "metro N=64 T=12 |phi(1)|^2 = " << mp.psi2_at_1 << " +/- " << mp.err_at_1 << "\n";
    std::cout << "<x^2> = " << mp.mean_x2 << " +/- " << mp.mean_x2_err << "   gap = " << gap << "\n";
    std::cout << "wrote " << dir.string() << " in " << sec << " s\n";
    return 0;
}

static void json_doubles(std::ostream& o, const std::vector<double>& v, int stride) {
    o << "[";
    bool first = true;
    for (size_t i = 0; i < v.size(); i += (size_t)std::max(1, stride)) {
        if (!first) o << ",";
        first = false;
        o << v[i];
    }
    if (!v.empty() && (v.size() - 1) % (size_t)std::max(1, stride) != 0) {
        if (!first) o << ",";
        o << v.back();
    }
    o << "]";
}

static int run_json(int argc, char** argv) {
    Potential pot = parse_potential(arg(argc, argv, "--pot", "ho"));
    MetroConfig c = base_cfg();
    c.pot = pot;
    c.N = argi(argc, argv, "--N", 36);
    c.T = argd(argc, argv, "--T", 8.0);
    c.mass = argd(argc, argv, "--mass", 1.0);
    c.hbar = argd(argc, argv, "--hbar", 1.0);
    c.sweeps = argi(argc, argv, "--sweeps", 8000);
    c.therm = argi(argc, argv, "--therm", 800);
    c.tune = argi(argc, argv, "--tune", 400);
    c.chains = argi(argc, argv, "--chains", 2);
    c.seed = (uint64_t)argi(argc, argv, "--seed", 1);
    c.bins = argi(argc, argv, "--bins", 140);
    c.xmin = argd(argc, argv, "--xmin", -6.0);
    c.xmax = argd(argc, argv, "--xmax", 6.0);
    c.batches = 8;
    c.save_paths = argi(argc, argv, "--paths", 8);
    if (c.N < 4) c.N = 4;
    if (c.N > 96) c.N = 96;
    if (c.sweeps < 400) c.sweeps = 400;
    if (c.sweeps > 40000) c.sweeps = 40000;
    if (c.sweeps < c.batches) c.batches = std::max(2, c.sweeps / 50);
    Spectrum sp = solve_schrodinger(pot, c.mass, c.hbar, c.xmin, c.xmax, 1201);
    MetroResult r = metropolis_ground_state(c);
    if (!r.ok || !sp.ok) {
        std::cout << "{\"ok\":false}\n";
        return 1;
    }
    int se_stride = std::max(1, (int)sp.x.size() / 360);
    std::cout << std::setprecision(8);
    std::cout << "{\"ok\":true,\"pot\":\"" << pot.name() << "\",\"label\":\"" << pot.label()
              << "\",\"N\":" << r.N << ",\"T\":" << r.T << ",\"eps\":" << r.eps
              << ",\"mass\":" << c.mass << ",\"hbar\":" << c.hbar
              << ",\"energy\":" << sp.energy
              << ",\"psi2_at_1\":" << r.psi2_at_1
              << ",\"psi2_at_1_err\":" << r.err_at_1
              << ",\"se_psi2_at_1\":" << interp(sp.x, sp.psi2, 1.0)
              << ",\"mean_x2\":" << r.mean_x2
              << ",\"acceptance\":" << r.acceptance
              << ",\"step\":" << r.step
              << ",\"x\":";
    json_doubles(std::cout, r.x, 1);
    std::cout << ",\"psi2\":";
    json_doubles(std::cout, r.psi2, 1);
    std::cout << ",\"err\":";
    json_doubles(std::cout, r.err, 1);
    std::cout << ",\"se_x\":";
    json_doubles(std::cout, sp.x, se_stride);
    std::cout << ",\"se_psi2\":";
    json_doubles(std::cout, sp.psi2, se_stride);
    std::cout << ",\"paths\":[";
    for (size_t p = 0; p < r.paths.size(); ++p) {
        if (p) std::cout << ",";
        json_doubles(std::cout, r.paths[p], 1);
    }
    std::cout << "]}\n";
    return 0;
}

int main(int argc, char** argv) {
    std::cout.setf(std::ios::fmtflags(0), std::ios::floatfield);
    std::cout.precision(10);
    if (argc < 2) {
        usage();
        return 1;
    }
    std::string cmd = argv[1];
    if (cmd == "selftest") return selftest();
    if (cmd == "suite") {
        fs::path dir = (argc >= 3) ? fs::path(argv[2]) : fs::path("data");
        return suite(dir);
    }
    if (cmd == "run") return run_json(argc, argv);
    usage();
    return 1;
}
