#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include "raylib.h"

#include "feynman/direct.hpp"
#include "feynman/exact_ho.hpp"
#include "feynman/metropolis.hpp"
#include "feynman/potential.hpp"
#include "feynman/histogram.hpp"
#include "feynman/research.hpp"
#include "feynman/schrodinger.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using namespace feynman;

static const int PAGE_MILE = 0;
static const int PAGE_CONV = 1;
static const int PAGE_PATHS = 2;
static const int PAGE_DIAG = 3;
static const int PAGE_GAP = 4;
static const int PAGE_CLASS = 5;
static const int PAGE_HBAR = 6;
static const int PAGE_POTS = 7;
static const int PAGE_DESIGN = 8;
static const int PAGE_RESEARCH = 9;
static const int PAGE_LIBRARY = 10;
static const int PAGE_COUNT = 11;

static const char* kPages[PAGE_COUNT] = {
    "Milestone", "Convergence", "Paths", "Diagnostics", "Gap", "Classical", "Planck", "Potentials", "Design", "Research", "Library"
};

static const char* kPotOrder[] = {"anharm", "doublewell", "morse", "quartic"};
static const char* kRunPots[] = {"ho", "anharm", "doublewell", "morse", "quartic"};
static const char* kRunLabels[] = {"harmonic", "anharmonic", "double well", "Morse", "quartic"};

static const Color COL_BG{16, 20, 28, 255};
static const Color COL_PANEL{26, 32, 44, 255};
static const Color COL_INSET{18, 22, 32, 255};
static const Color COL_INK{243, 239, 230, 255};
static const Color COL_MUTED{158, 152, 140, 255};
static const Color COL_LINE{54, 64, 82, 255};
static const Color COL_TEAL{32, 148, 140, 255};
static const Color COL_TEAL_DIM{32, 148, 140, 70};
static const Color COL_COPPER{204, 92, 58, 255};
static const Color COL_GOLD{214, 168, 72, 255};
static const Color COL_NAVY{120, 156, 196, 255};
static const Color COL_HOT{42, 52, 70, 255};

static Font g_ui{};
static Font g_code{};

static Font font_or_default(Font font) {
    if (font.texture.id == 0) return GetFontDefault();
    return font;
}

static void text(Font font, const char* s, float x, float y, float size, Color c) {
    DrawTextEx(font_or_default(font), s, Vector2{x, y}, size, 0.4f, c);
}

static float text_width(Font font, const char* s, float size) {
    return MeasureTextEx(font_or_default(font), s, size, 0.4f).x;
}

struct Pack {
    int generation = 0;
    int N = 7;
    double T = 4.0;
    double mass = 1.0;
    double hbar_run = 1.0;
    const char* run_pot = "ho";
    DirectResult direct;
    Spectrum schrodinger;
    HoLattice lattice;
    MetroResult metro_small;
    MetroResult metro;
    std::vector<BridgePath> bridges;
    int Ns[7] = {4, 7, 12, 24, 48, 96, 160};
    double Ts[6] = {2, 3, 4, 8, 12, 16};
    double psi1[7][6]{};
    double alpha[7][6]{};
    double hbars[6] = {0.4, 0.7, 1.0, 1.6, 2.4, 3.5};
    MetroResult hbar[6];
    const char* pot_name[4] = {"anharm", "doublewell", "morse", "quartic"};
    MetroResult pot_metro[4];
    Spectrum pot_se[4];
    std::vector<double> acf;
    double gap = 0;
    double forbid = 0;
    double analytic1 = 0;
    double ho_e0 = 0;
    double ho_e1 = 0;
    double ho_energy = 0;
    double ho_energy_err = 0;
    double ho_cross = 0;
    double dw_e0 = 0;
    double dw_e1 = 0;
    double dw_energy = 0;
    double dw_cross = 0;
    double split_h[7] = {0.45, 0.60, 0.80, 1.00, 1.30, 1.70, 2.20};
    double split_gap[7]{};
    double res_h[8] = {0.35, 0.50, 0.70, 1.00, 1.40, 2.00, 2.80, 3.60};
    double res_exact[8]{};
    double res_mc[6]{};
    double res_bar[7]{};
    double an_e0 = 0, an_gap = 0;
    double mo_e0 = 0, mo_gap = 0;
    double qu_e0 = 0, qu_gap = 0;
    double ms_se = 0, ms_direct = 0, ms_metro = 0;
    int direct_N = 0;
    LargeHist hist;
};

static std::mutex g_mu;
static std::condition_variable g_cv;
static std::shared_ptr<Pack> g_pack;
static int g_job_seq = 0;
static int g_job_kind = 1;
static int g_job_N = 7;
static double g_job_T = 4;
static double g_job_mass = 1;
static double g_job_hbar = 1;
static int g_job_pot = 0;
static bool g_quit = false;
static double g_set_N = 7;
static double g_set_T = 4;
static double g_set_hbar = 1;
static double g_set_mass = 1;
static int g_set_pot = 0;
static std::atomic<int> g_stage{0};
static std::atomic<bool> g_busy{true};

static const char* stage_name(int s) {
    switch (s) {
    case 1: return "Solving the Schrodinger equation";
    case 2: return "Summing random paths on the milestone lattice";
    case 3: return "Exact lattice integral, every N and T";
    case 4: return "Metropolis on the same small lattice";
    case 5: return "Long Metropolis chain";
    case 6: return "Scanning hbar";
    case 7: return "Anharmonic, double well, Morse, quartic";
    case 8: return "Drawing the pinned bridges";
    case 10: return "Histogram of the slices";
    default: return "Starting the engine";
    }
}

static MetroConfig mc_base() {
    MetroConfig c{};
    c.mass = 1.0;
    c.hbar = 1.0;
    c.pot = parse_potential("ho");
    c.step = 0.0;
    c.x0 = 0.0;
    c.batches = 8;
    c.do_corr = false;
    c.do_world = false;
    c.do_trace = false;
    c.save_paths = 0;
    c.seed = 1;
    c.tune = 300;
    c.therm = 500;
    c.sweeps = 6000;
    c.chains = 2;
    c.bins = 140;
    c.xmin = -6;
    c.xmax = 6;
    c.N = 32;
    c.T = 8;
    return c;
}

static void pot_window(const char* name, double& xmin, double& xmax) {
    if (strcmp(name, "doublewell") == 0) { xmin = -3.4; xmax = 3.4; }
    else if (strcmp(name, "morse") == 0) { xmin = -1.6; xmax = 6.0; }
    else if (strcmp(name, "anharm") == 0) { xmin = -3.6; xmax = 3.6; }
    else if (strcmp(name, "quartic") == 0) { xmin = -2.6; xmax = 2.6; }
    else { xmin = -5.5; xmax = 5.5; }
}

static void fill_milestone(Pack& p, int N, double T, double mass, double hbar, int poti) {
    if (N < 4) N = 4;
    if (N > 48) N = 48;
    if (T < 0.5) T = 0.5;
    if (!(mass > 0.0)) mass = 1.0;
    if (!(hbar > 0.0)) hbar = 1.0;
    if (poti < 0 || poti > 4) poti = 0;
    const char* pname = kRunPots[poti];
    p.N = N;
    p.T = T;
    p.mass = mass;
    p.hbar_run = hbar;
    p.run_pot = pname;
    Potential pot = parse_potential(pname);
    bool ho = poti == 0;
    double xmin = -5.5, xmax = 5.5;
    pot_window(pname, xmin, xmax);
    g_stage = 1;
    p.schrodinger = solve_schrodinger(pot, mass, hbar, ho ? -10.0 : xmin, ho ? 10.0 : xmax, 1401);
    p.analytic1 = ho ? ho_ground_psi2(1.0, mass, hbar, 1.0) : 0.0;
    g_stage = 2;
    DirectConfig dc{};
    dc.N = N;
    dc.T = T;
    dc.mass = mass;
    dc.hbar = hbar;
    dc.pot = pot;
    dc.xmin = xmin;
    dc.xmax = xmax;
    dc.dx = 0.05;
    dc.samples = (N <= 12) ? 40000 : 20000;
    dc.batches = 8;
    dc.seed = 42;
    dc.threads = 4;
    p.direct = direct_ground_density(dc);
    p.direct_N = N;
    if (ho) p.lattice = lattice_alpha(N, T, mass, hbar, 1.0);
    else p.lattice = HoLattice{};
    g_stage = 4;
    MetroConfig small = mc_base();
    small.pot = pot;
    small.N = N;
    small.T = T;
    small.mass = mass;
    small.hbar = hbar;
    small.sweeps = 5000;
    small.chains = 2;
    small.batches = 6;
    small.bins = 140;
    small.xmin = xmin;
    small.xmax = xmax;
    small.seed = 3;
    p.metro_small = metropolis_ground_state(small);
    g_stage = 8;
    int bn = std::min(N, 16);
    p.bridges = sample_bridges(1.0, bn, T * (double)bn / (double)N, mass, hbar, pot, 20, 99);
}

static std::vector<double> acf_of(const std::vector<double>& series, const std::vector<int>& flag, int lags);
static double fit_gap(const MetroResult& R);
static double tail_outside(const MetroResult& R);
static bool job_stale(int seen);
static void fill_quick(Pack& p, int N, double T, double mass, double hbar, int poti);
static void fill_large(Pack& p, int N, double T, double mass, double hbar, int poti);

static void fill_live_chain(Pack& p) {
    g_stage = 5;
    const char* pname = p.run_pot ? p.run_pot : "ho";
    Potential pot = parse_potential(pname);
    double xmin = -6, xmax = 6;
    pot_window(pname, xmin, xmax);
    MetroConfig prod = mc_base();
    prod.pot = pot;
    prod.N = p.N;
    prod.T = p.T;
    prod.mass = p.mass;
    prod.hbar = p.hbar_run;
    prod.sweeps = 5000;
    prod.therm = 500;
    prod.tune = 300;
    prod.chains = 2;
    prod.batches = 6;
    prod.bins = 140;
    prod.xmin = xmin;
    prod.xmax = xmax;
    prod.seed = 5;
    prod.save_paths = 8;
    prod.do_corr = true;
    prod.do_world = true;
    prod.do_trace = true;
    p.metro = metropolis_ground_state(prod);
    p.acf = acf_of(p.metro.trace_x2, p.metro.trace_measuring, 60);
    p.gap = fit_gap(p.metro);
    p.forbid = tail_outside(p.metro);
}

static std::vector<double> acf_of(const std::vector<double>& series, const std::vector<int>& flag, int lags) {
    std::vector<double> x;
    for (size_t i = 0; i < series.size(); ++i) {
        if (i < flag.size() && flag[i] == 0) continue;
        x.push_back(series[i]);
    }
    std::vector<double> out(lags, 0.0);
    if (x.size() < 8) return out;
    double mean = 0.0;
    for (double v : x) mean += v;
    mean /= (double)x.size();
    double var = 0.0;
    for (double v : x) var += (v - mean) * (v - mean);
    var /= (double)x.size();
    if (!(var > 1e-12)) return out;
    for (int lag = 0; lag < lags; ++lag) {
        double s = 0.0;
        int n = (int)x.size() - lag;
        for (int i = 0; i < n; ++i) s += (x[i] - mean) * (x[i + lag] - mean);
        out[lag] = (s / n) / var;
    }
    return out;
}

static double fit_gap(const MetroResult& R) {
    double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 0; i < R.corr_tau.size(); ++i) {
        double tau = R.corr_tau[i];
        double c = R.corr_xx[i];
        if (tau < 0.35 || tau > 2.4 || !(c > 0.0)) continue;
        double y = std::log(c);
        sx += tau;
        sy += y;
        sxx += tau * tau;
        sxy += tau * y;
        n += 1;
    }
    double den = n * sxx - sx * sx;
    if (n < 3 || !(std::fabs(den) > 0)) return 0;
    return -(n * sxy - sx * sy) / den;
}

static double tail_outside(const MetroResult& R) {
    if (R.x.size() < 2) return 0;
    double p = 0;
    for (size_t i = 0; i + 1 < R.x.size(); ++i) {
        double mid = 0.5 * (R.x[i] + R.x[i + 1]);
        if (std::fabs(mid) > 1.0) p += 0.5 * (R.psi2[i] + R.psi2[i + 1]) * (R.x[i + 1] - R.x[i]);
    }
    return p;
}

static void fill_gallery(Pack& p) {
    Potential ho = parse_potential("ho");
    g_stage = 3;
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 6; ++j) {
            HoLattice a = lattice_alpha(p.Ns[i], p.Ts[j], 1.0, 1.0, 1.0);
            p.alpha[i][j] = a.ok ? a.alpha : 0;
            p.psi1[i][j] = a.ok ? gaussian_psi2(1.0, a.alpha) : 0;
        }
    }
    g_stage = 5;
    MetroConfig prod = mc_base();
    prod.N = 48;
    prod.T = 12.0;
    prod.sweeps = 8000;
    prod.therm = 800;
    prod.tune = 400;
    prod.chains = 3;
    prod.batches = 8;
    prod.bins = 160;
    prod.xmin = -6;
    prod.xmax = 6;
    prod.seed = 5;
    prod.save_paths = 10;
    prod.do_corr = true;
    prod.do_world = true;
    prod.do_trace = true;
    p.metro = metropolis_ground_state(prod);
    p.acf = acf_of(p.metro.trace_x2, p.metro.trace_measuring, 60);
    p.gap = fit_gap(p.metro);
    p.forbid = tail_outside(p.metro);
    g_stage = 6;
    for (int i = 0; i < 6; ++i) {
        MetroConfig c = mc_base();
        c.N = 32;
        c.T = 8.0;
        c.hbar = p.hbars[i];
        c.sweeps = 4000;
        c.chains = 2;
        c.bins = 120;
        c.xmin = -7;
        c.xmax = 7;
        c.seed = 100 + i * 17;
        p.hbar[i] = metropolis_ground_state(c);
    }
    g_stage = 7;
    double xmin[4] = {-3.4, -3.2, -1.5, -2.4};
    double xmax[4] = {3.4, 3.2, 5.6, 2.4};
    for (int i = 0; i < 4; ++i) {
        Potential pot = parse_potential(kPotOrder[i]);
        p.pot_name[i] = kPotOrder[i];
        p.pot_se[i] = solve_schrodinger(pot, 1.0, 1.0, xmin[i], xmax[i], 1201);
        MetroConfig c = mc_base();
        c.pot = pot;
        c.N = 36;
        c.T = 8.0;
        c.sweeps = 5000;
        c.chains = 2;
        c.bins = 120;
        c.xmin = xmin[i];
        c.xmax = xmax[i];
        c.seed = 200 + i * 13;
        c.do_corr = (i == 1);
        p.pot_metro[i] = metropolis_ground_state(c);
    }
    p.ho_e0 = p.schrodinger.energy;
    Spectrum ho1 = solve_schrodinger(ho, 1.0, 1.0, -10.0, 10.0, 1401, 1);
    p.ho_e1 = ho1.energy;
    p.ho_energy = p.metro.mean_energy;
    p.ho_energy_err = p.metro.mean_energy_err;
    p.ho_cross = p.metro.mean_crossings;
    p.dw_e0 = p.pot_se[1].energy;
    Potential dw = parse_potential("doublewell");
    Spectrum dw1 = solve_schrodinger(dw, 1.0, 1.0, -3.2, 3.2, 1201, 1);
    p.dw_e1 = dw1.energy;
    p.dw_energy = p.pot_metro[1].mean_energy;
    p.dw_cross = p.pot_metro[1].mean_crossings;
    for (int i = 0; i < 7; ++i) {
        Spectrum a = solve_schrodinger(dw, 1.0, p.split_h[i], -3.4, 3.4, 1001, 0);
        Spectrum b = solve_schrodinger(dw, 1.0, p.split_h[i], -3.4, 3.4, 1001, 1);
        p.split_gap[i] = (a.ok && b.ok) ? (b.energy - a.energy) : 0.0;
    }
    for (int i = 0; i < 8; ++i) p.res_exact[i] = ho_outside_fixed(p.res_h[i], 1.0);
    for (int i = 0; i < 6; ++i) p.res_mc[i] = tail_outside(p.hbar[i]);
    Potential dwpot = parse_potential("doublewell");
    for (int i = 0; i < 7; ++i) {
        Spectrum g = solve_schrodinger(dwpot, 1.0, p.split_h[i], -3.4, 3.4, 801, 0);
        double xin = 0, xout = 0;
        if (g.ok && double_well_edges(g.energy, dwpot.a, dwpot.b, xin, xout)) p.res_bar[i] = density_between(g, -xin, xin);
    }
    auto gap_of = [](const char* name, double xmin, double xmax, double& e0, double& gap) {
        Potential pot = parse_potential(name);
        Spectrum a = solve_schrodinger(pot, 1.0, 1.0, xmin, xmax, 1201, 0);
        Spectrum b = solve_schrodinger(pot, 1.0, 1.0, xmin, xmax, 1201, 1);
        e0 = a.energy;
        gap = (a.ok && b.ok) ? b.energy - a.energy : 0.0;
    };
    gap_of("anharm", -3.6, 3.6, p.an_e0, p.an_gap);
    gap_of("morse", -1.6, 6.0, p.mo_e0, p.mo_gap);
    gap_of("quartic", -2.6, 2.6, p.qu_e0, p.qu_gap);
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    solve_schrodinger(ho, 1.0, 1.0, -10.0, 10.0, 801, 0);
    auto t1 = clock::now();
    DirectConfig dc{};
    dc.N = 12;
    dc.T = 4;
    dc.mass = 1;
    dc.hbar = 1;
    dc.pot = ho;
    dc.xmin = 0;
    dc.xmax = 2;
    dc.dx = 0.1;
    dc.samples = 2000;
    dc.batches = 4;
    dc.seed = 3;
    dc.threads = 4;
    direct_ground_density(dc);
    auto t2 = clock::now();
    MetroConfig timed = mc_base();
    timed.N = 16;
    timed.T = 4;
    timed.sweeps = 800;
    timed.therm = 100;
    timed.tune = 80;
    timed.chains = 1;
    timed.batches = 4;
    timed.bins = 60;
    timed.xmin = -4;
    timed.xmax = 4;
    timed.do_corr = false;
    timed.do_world = false;
    timed.do_trace = false;
    metropolis_ground_state(timed);
    auto t3 = clock::now();
    auto millis = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    p.ms_se = millis(t0, t1);
    p.ms_direct = millis(t1, t2);
    p.ms_metro = millis(t2, t3);
}

static void worker_main() {
    int seen = 0;
    while (true) {
        int kind, N, poti;
        double T, mass, hbar;
        {
            std::unique_lock<std::mutex> lock(g_mu);
            g_cv.wait(lock, [&] { return g_quit || g_job_seq != seen; });
            if (g_quit) return;
            kind = g_job_kind;
            N = g_job_N;
            T = g_job_T;
            mass = g_job_mass;
            hbar = g_job_hbar;
            poti = g_job_pot;
            seen = g_job_seq;
        }
        g_busy = true;
        auto next = std::make_shared<Pack>();
        if (kind != 1) {
            std::shared_ptr<Pack> old;
            {
                std::lock_guard<std::mutex> lock(g_mu);
                old = g_pack;
            }
            if (old) *next = *old;
        }
        bool publish = true;
        if (kind == 1) {
            fill_milestone(*next, N, T, mass, hbar, poti);
            fill_gallery(*next);
        } else if (kind == 2) {
            fill_milestone(*next, N, T, mass, hbar, poti);
            fill_live_chain(*next);
        } else if (kind == 4) {
            if (!job_stale(seen)) {
                fill_large(*next, N, T, mass, hbar, poti);
                if (job_stale(seen)) publish = false;
            } else {
                publish = false;
            }
        } else if (!job_stale(seen)) {
            fill_quick(*next, N, T, mass, hbar, poti);
            if (job_stale(seen)) publish = false;
        } else {
            publish = false;
        }
        if (publish) {
            std::lock_guard<std::mutex> lock(g_mu);
            if ((kind == 3 || kind == 4) && g_job_seq != seen) publish = false;
            else {
                next->generation = g_pack ? g_pack->generation + 1 : 1;
                g_pack = next;
            }
        }
        g_stage = 9;
        g_busy = false;
    }
}

static void request_job(int kind, int N, double T, double mass, double hbar, int poti) {
    {
        std::lock_guard<std::mutex> lock(g_mu);
        g_job_kind = kind;
        g_job_N = N;
        g_job_T = T;
        g_job_mass = mass;
        g_job_hbar = hbar;
        g_job_pot = poti;
        g_job_seq += 1;
    }
    g_cv.notify_one();
}

static void launch_run() {
    if (g_busy) return;
    int N = (int)std::lround(g_set_N);
    if (N < 4) N = 4;
    if (N > 1000000) N = 1000000;
    int kind = (N > 48) ? 4 : 2;
    request_job(kind, N, g_set_T, g_set_mass, g_set_hbar, g_set_pot);
}

struct Live {
    int N = 7;
    double T = 4;
    double mass = 1;
    double hbar = 1;
    int pot = 0;
    Spectrum se;
    HoLattice lattice{};
    bool ready = false;
};

static Live g_live;
static bool g_play = false;

static bool job_stale(int seen) {
    std::lock_guard<std::mutex> lock(g_mu);
    return g_quit || g_job_seq != seen;
}

static void refresh_exact() {
    int N = (int)std::lround(g_set_N);
    if (N < 4) N = 4;
    if (N > 1000000) N = 1000000;
    double T = g_set_T;
    double mass = g_set_mass;
    double hbar = g_set_hbar;
    int poti = g_set_pot;
    if (poti < 0 || poti > 4) poti = 0;
    if (T < 0.5) T = 0.5;
    if (!(mass > 0.0)) mass = 1.0;
    if (!(hbar > 0.0)) hbar = 1.0;
    const char* pname = kRunPots[poti];
    Potential pot = parse_potential(pname);
    bool ho = poti == 0;
    double xmin = -5.5, xmax = 5.5;
    pot_window(pname, xmin, xmax);
    g_live.N = N;
    g_live.T = T;
    g_live.mass = mass;
    g_live.hbar = hbar;
    g_live.pot = poti;
    g_live.se = solve_schrodinger(pot, mass, hbar, ho ? -10.0 : xmin, ho ? 10.0 : xmax, 801);
    if (ho && N <= 48) g_live.lattice = lattice_alpha(N, T, mass, hbar, 1.0);
    else g_live.lattice = HoLattice{};
    g_live.ready = true;
}

static void fill_quick(Pack& p, int N, double T, double mass, double hbar, int poti) {
    if (N < 4) N = 4;
    if (N > 48) N = 48;
    if (T < 0.5) T = 0.5;
    if (!(mass > 0.0)) mass = 1.0;
    if (!(hbar > 0.0)) hbar = 1.0;
    if (poti < 0 || poti > 4) poti = 0;
    const char* pname = kRunPots[poti];
    p.N = N;
    p.T = T;
    p.mass = mass;
    p.hbar_run = hbar;
    p.run_pot = pname;
    Potential pot = parse_potential(pname);
    bool ho = poti == 0;
    double xmin = 0.0, xmax = 2.0;
    if (!ho) pot_window(pname, xmin, xmax);
    p.analytic1 = ho ? ho_ground_psi2(1.0, mass, hbar, 1.0) : 0.0;
    p.schrodinger = solve_schrodinger(pot, mass, hbar, ho ? -10.0 : xmin, ho ? 10.0 : xmax, 801);
    if (ho) p.lattice = lattice_alpha(N, T, mass, hbar, 1.0);
    else p.lattice = HoLattice{};
    DirectConfig dc{};
    dc.N = N;
    dc.T = T;
    dc.mass = mass;
    dc.hbar = hbar;
    dc.pot = pot;
    dc.xmin = xmin;
    dc.xmax = xmax;
    dc.dx = ho ? 0.05 : 0.1;
    dc.samples = 2500;
    dc.batches = 4;
    dc.seed = 42 + (uint64_t)N * 17u;
    dc.threads = 4;
    p.direct = direct_ground_density(dc);
    p.direct_N = N;
}

static void fill_large(Pack& p, int N, double T, double mass, double hbar, int poti) {
    if (N < 4) N = 4;
    if (N > 1000000) N = 1000000;
    if (T < 0.5) T = 0.5;
    if (!(mass > 0.0)) mass = 1.0;
    if (!(hbar > 0.0)) hbar = 1.0;
    if (poti < 0 || poti > 4) poti = 0;
    const char* pname = kRunPots[poti];
    g_stage = 10;
    Potential pot = parse_potential(pname);
    double xmin = -6.0, xmax = 6.0;
    if (poti != 0) pot_window(pname, xmin, xmax);
    p.hist = large_histogram(pot, N, T, mass, hbar, 400, 0, 100, xmin, xmax, 11u + (uint64_t)N);
    p.N = N;
    p.T = T;
    p.mass = mass;
    p.hbar_run = hbar;
    p.run_pot = pname;
    if (p.hist.ok) write_histogram_csv("data/histogram.csv", p.hist, pname, mass, hbar);
}

static void kick_live() {
    int N = (int)std::lround(g_set_N);
    if (N > 48) return;
    request_job(3, N, g_set_T, g_set_mass, g_set_hbar, g_set_pot);
}

static bool sample_matches(const Pack& p) {
    if (!g_live.ready || !p.direct.ok) return false;
    if (p.direct_N != g_live.N) return false;
    if (p.N != g_live.N) return false;
    if (std::fabs(p.T - g_live.T) > 1e-9) return false;
    if (std::fabs(p.mass - g_live.mass) > 1e-9) return false;
    if (std::fabs(p.hbar_run - g_live.hbar) > 1e-9) return false;
    if (!p.run_pot || strcmp(p.run_pot, kRunPots[g_live.pot]) != 0) return false;
    return true;
}

static bool hist_matches(const Pack& p) {
    if (!g_live.ready || !p.hist.ok || g_live.N <= 48) return false;
    if (p.hist.N != g_live.N) return false;
    if (std::fabs(p.hist.T - g_live.T) > 1e-6) return false;
    if (std::fabs(p.mass - g_live.mass) > 1e-9) return false;
    if (std::fabs(p.hbar_run - g_live.hbar) > 1e-9) return false;
    if (!p.run_pot || strcmp(p.run_pot, kRunPots[g_live.pot]) != 0) return false;
    return true;
}

struct Axis {
    Rectangle box;
    double x0, x1, y0, y1;
    float L, R, T, B;
    Vector2 map(double x, double y) const {
        float px = box.x + L + (float)((x - x0) / (x1 - x0)) * (box.width - L - R);
        float py = box.y + box.height - B - (float)((y - y0) / (y1 - y0)) * (box.height - T - B);
        return Vector2{px, py};
    }
    Rectangle plot_rect() const {
        return Rectangle{box.x + L, box.y + T, box.width - L - R, box.height - T - B};
    }
};

static Axis make_axis(Rectangle box, double x0, double x1, double y0, double y1) {
    Axis a;
    a.box = box;
    a.x0 = x0;
    a.x1 = x1;
    a.y0 = y0;
    a.y1 = y1;
    a.L = 54;
    a.R = 16;
    a.T = 16;
    a.B = 36;
    return a;
}

static void draw_frame(const Axis& a, const char* xlab, const char* ylab) {
    Rectangle inner = a.plot_rect();
    DrawRectangleRec(inner, COL_INSET);
    DrawRectangleLinesEx(inner, 1.0f, COL_LINE);
    for (int i = 1; i <= 4; ++i) {
        float y = inner.y + inner.height * (float)i / 5.0f;
        DrawLineEx(Vector2{inner.x, y}, Vector2{inner.x + inner.width, y}, 1.0f, Color{255, 255, 255, 16});
    }
    char buf[32];
    for (int i = 0; i <= 4; ++i) {
        double x = a.x0 + (a.x1 - a.x0) * i / 4.0;
        if (std::fabs(x) < 1e-6) x = 0.0;
        snprintf(buf, sizeof buf, "%.4g", x);
        Vector2 p = a.map(x, a.y0);
        float w = text_width(g_ui, buf, 14);
        text(g_ui, buf, p.x - w * 0.5f, a.box.y + a.box.height - 22, 14, COL_MUTED);
    }
    for (int i = 0; i <= 3; ++i) {
        double y = a.y0 + (a.y1 - a.y0) * i / 3.0;
        snprintf(buf, sizeof buf, "%.2f", y);
        Vector2 p = a.map(a.x0, y);
        text(g_ui, buf, a.box.x + 4, p.y - 8, 13, COL_MUTED);
    }
    float lw = text_width(g_ui, xlab, 15);
    text(g_ui, xlab, inner.x + inner.width * 0.5f - lw * 0.5f, a.box.y + a.box.height - 4, 15, COL_MUTED);
    (void)ylab;
}

static void stroke(const Axis& a, const std::vector<double>& xs, const std::vector<double>& ys, Color c, float thick, bool dashed) {
    if (xs.size() < 2 || ys.size() < 2) return;
    BeginScissorMode((int)a.plot_rect().x, (int)a.plot_rect().y, (int)a.plot_rect().width, (int)a.plot_rect().height);
    int run = 0;
    Vector2 prev{};
    bool has = false;
    size_t n = std::min(xs.size(), ys.size());
    for (size_t i = 0; i < n; ++i) {
        if (xs[i] < a.x0 || xs[i] > a.x1) {
            has = false;
            continue;
        }
        Vector2 p = a.map(xs[i], ys[i]);
        if (has) {
            if (!dashed || ((run / 4) % 2 == 0)) DrawLineEx(prev, p, thick, c);
            run++;
        }
        prev = p;
        has = true;
    }
    EndScissorMode();
}

static void band(const Axis& a, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& err, Color c) {
    if (xs.size() < 2) return;
    BeginScissorMode((int)a.plot_rect().x, (int)a.plot_rect().y, (int)a.plot_rect().width, (int)a.plot_rect().height);
    Vector2 prev_lo{}, prev_hi{};
    bool has = false;
    size_t n = std::min(xs.size(), std::min(ys.size(), err.size()));
    for (size_t i = 0; i < n; ++i) {
        if (xs[i] < a.x0 || xs[i] > a.x1) {
            has = false;
            continue;
        }
        Vector2 lo = a.map(xs[i], ys[i] - err[i]);
        Vector2 hi = a.map(xs[i], ys[i] + err[i]);
        if (has) {
            DrawTriangle(prev_lo, prev_hi, hi, c);
            DrawTriangle(prev_lo, hi, lo, c);
        }
        prev_lo = lo;
        prev_hi = hi;
        has = true;
    }
    EndScissorMode();
}

static void legend(Rectangle host, const char** names, const Color* colors, int n) {
    float y = host.y + 8;
    float x = host.x + host.width - 210;
    DrawRectangleRounded(Rectangle{x - 8, y - 4, 200, 18.0f * n + 8}, 0.15f, 4, Color{16, 20, 28, 210});
    for (int i = 0; i < n; ++i) {
        DrawRectangleRec(Rectangle{x, y + 4, 14, 4}, colors[i]);
        text(g_ui, names[i], x + 20, y, 14, COL_INK);
        y += 18;
    }
}

static bool button(Rectangle r, const char* label, bool on) {
    Vector2 m = GetMousePosition();
    bool hot = CheckCollisionPointRec(m, r);
    Color bg = on ? COL_TEAL : (hot ? COL_HOT : Color{32, 40, 54, 255});
    DrawRectangleRounded(r, 0.4f, 6, bg);
    if (on) DrawRectangleRoundedLines(r, 0.4f, 6, COL_TEAL);
    float w = text_width(g_ui, label, 16);
    text(g_ui, label, r.x + (r.width - w) * 0.5f, r.y + (r.height - 16) * 0.5f, 16, on ? Color{12, 24, 24, 255} : COL_INK);
    return hot && IsMouseButtonPressed(MOUSE_BUTTON_LEFT);
}

static void card(Rectangle r) {
    DrawRectangleRounded(r, 0.04f, 8, COL_PANEL);
}

static const char* kScriptMile =
    "direct_ground_density          direct.hpp\n"
    "\n"
    "eps = T / N\n"
    "for x in the grid from 0 to 2:\n"
    "    weights = []\n"
    "    repeat samples:\n"
    "        path = brownian_bridge(x -> x)\n"
    "        S = sum eps * V(midpoint)\n"
    "        weights.append(exp(-S))\n"
    "    G(x) = mean(weights)\n"
    "|phi|^2 = G / integral(G)\n"
    "report |phi(1)|^2 and the batch error\n"
    "\n"
    "The bridge samples the kinetic piece\n"
    "exactly. The weight is only the\n"
    "potential. Together that is e^(-eps E)\n"
    "from the brief, equation 14.8.\n"
    "\n"
    "link_action(left, right) is the one\n"
    "term inside the sum for E.";

static const char* kScriptConv =
    "lattice_alpha                  exact_ho.hpp\n"
    "\n"
    "For the oscillator the integral is\n"
    "Gaussian, so it can be done with a\n"
    "Cholesky factor and no random numbers.\n"
    "\n"
    "for N in 4, 7, 12, 24, 48, 96, 160:\n"
    "    for T in 2, 3, 4, 8, 12, 16:\n"
    "        alpha = lattice_alpha(N, T)\n"
    "        |phi(1)|^2 = sqrt(alpha/pi)\n"
    "                     * exp(-alpha)\n"
    "\n"
    "T = 3 is the other time written in\n"
    "the brief. T = 4 is the one next to\n"
    "eps = T/N. Both are on this grid.\n"
    "\n"
    "Short T is the error that matters.\n"
    "Adding slices at fixed T = 4 does\n"
    "not reach the ground state.";

static const char* kScriptPaths =
    "metropolis_ground_state     metropolis.hpp\n"
    "\n"
    "path = N positions on a ring\n"
    "for sweep in sweeps:\n"
    "    for each site k:\n"
    "        propose x' = x[k] + random(-d, d)\n"
    "        dS = new link_action - old\n"
    "             on the two neighbouring links\n"
    "        keep x' if dS <= 0\n"
    "        else keep it if uniform() < exp(-dS)\n"
    "histogram of x  ->  |phi(x)|^2\n"
    "\n"
    "d is tuned on a pilot, then frozen.\n"
    "Four ideas from the brief are here:\n"
    "one chain versus several, the\n"
    "autocorrelation, the correlator\n"
    "<x(tau) x(0)>, and the paths themselves.";

static const char* kScriptDiag =
    "diagnostics from the long chain\n"
    "\n"
    "trace of <x^2> and the acceptance\n"
    "after the step size is frozen\n"
    "\n"
    "correlator\n"
    "  C(tau) = <x(tau) x(0)>\n"
    "  theory = 0.5 * cosh(T/2 - tau)\n"
    "                  / sinh(T/2)\n"
    "\n"
    "gap estimate: slope of log C(tau)\n"
    "between tau = 0.35 and 2.4\n"
    "The true gap of this oscillator is 1.\n"
    "The fit is rough because the periodic\n"
    "cosh is not a pure exponential.\n"
    "\n"
    "autocorrelation of <x^2> across sweeps\n"
    "is why the error uses batch means,\n"
    "not every sweep as if it were new.";

static const char* kScriptClass =
    "classical dwell time\n"
    "\n"
    "E = 1/2,  V = x^2 / 2\n"
    "turning points at x = +/- 1\n"
    "rho_cl(x) = 1 / (pi * sqrt(1 - x^2))\n"
    "for |x| < 1, and 0 outside.\n"
    "\n"
    "The quantum density is the long\n"
    "Metropolis histogram. The fraction\n"
    "past |x| = 1 is the tunnelling tail.\n"
    "A classical particle with this energy\n"
    "never gets there.\n"
    "\n"
    "hbar enters only as a divisor on the\n"
    "action. Turn it down on the Planck page.";

static const char* kScriptHbar =
    "same metropolis_ground_state\n"
    "only hbar changes\n"
    "\n"
    "S = sum  mass*(dx)^2 / (2*eps*hbar)\n"
    "       + eps * V(mid) / hbar\n"
    "\n"
    "for hbar in 0.4, 0.7, 1, 1.6, 2.4, 3.5:\n"
    "    run the chain\n"
    "    plot the histogram\n"
    "    overlay the exact ground state\n"
    "    sqrt(1/(pi*hbar)) * exp(-x^2/hbar)\n"
    "\n"
    "Small hbar squeezes the particle\n"
    "onto the origin. That is the classical\n"
    "limit the brief suggests studying.";

static const char* kScriptPots =
    "Potential::operator()         potential.hpp\n"
    "\n"
    "anharm     x^2/2 + x^4/4\n"
    "doublewell 0.8 * (x^2 - 2.25)^2\n"
    "morse      3 * (1 - exp(-0.65 x))^2\n"
    "quartic    0.25 * x^4\n"
    "\n"
    "The chain does not know which one.\n"
    "It only calls V(midpoint).\n"
    "\n"
    "To add another:\n"
    "  1. a case in potential.hpp\n"
    "  2. one line in kPotOrder in studio.cpp\n"
    "The Schrodinger solver and the\n"
    "Metropolis update stay as they are.\n"
    "Quartic is the extra one, past the\n"
    "three in the write-up.";

static const char* kScriptGap =
    "solve_schrodinger(..., level)   schrodinger.hpp\n"
    "energy_estimator                action.hpp\n"
    "sign_crossings\n"
    "\n"
    "level 0 is the ground state.\n"
    "level 1 is the first excited state.\n"
    "The gap is E1 - E0.\n"
    "\n"
    "On a path the same gap is the\n"
    "slope of <x(tau) x(0)>.\n"
    "That works for the oscillator.\n"
    "\n"
    "On the double well the splitting\n"
    "is small, and a local move almost\n"
    "never crosses the barrier.\n"
    "The correlator stays flat.\n"
    "The spectrum still has the gap.\n"
    "\n"
    "The curve is that gap against hbar.\n"
    "It falls off very fast as hbar drops.\n"
    "That is the tunnelling the year\n"
    "of work is about.";

static const char* kScriptResearch =
    "ho_outside_fixed            research.hpp\n"
    "ho_outside_turning\n"
    "double_well_edges\n"
    "density_between\n"
    "solve_schrodinger(..., level)\n"
    "\n"
    "Outside |x| = 1 the oscillator\n"
    "tail is erfc(1 / sqrt(hbar)).\n"
    "It dies as hbar drops.\n"
    "\n"
    "Outside the turning point the\n"
    "same tail is erfc(1), for every\n"
    "hbar. The wall moves out with\n"
    "the Gaussian.\n"
    "\n"
    "The double-well gap is not the\n"
    "amount of density in the barrier.\n"
    "At small hbar the barrier still\n"
    "holds about a tenth of the\n"
    "particle, and the gap is tiny.\n"
    "\n"
    "Anharmonic, Morse and quartic\n"
    "have no closed form. The gap is\n"
    "just the second eigenvalue.";

static const char* kScriptLibrary =
    "#include \"feynman/feynman.hpp\"\n"
    "\n"
    "link_action(...)\n"
    "direct_ground_density(...)\n"
    "metropolis_ground_state(...)\n"
    "solve_schrodinger(...)\n"
    "lattice_alpha(...)\n"
    "ho_outside_fixed(...)\n"
    "\n"
    "WRITEUP.md   the milestone\n"
    "THEORY.md    the tunnel splitting\n"
    "README.md    the build\n"
    "\n"
    "A new simulation is a new file\n"
    "that includes the header.\n"
    "It should not copy studio.cpp.\n"
    "\n"
    "Use Metropolis for a whole curve.\n"
    "The direct sum rebuilds a bridge\n"
    "at every x.\n"
    "\n"
    "The slow case is the double well.\n"
    "Local steps do not cross, so\n"
    "more sweeps do not help.\n"
    "The kink update is the fix.";

static const char* kScriptDesign =
    "DESIGN  (what the brief asks first)\n"
    "\n"
    "data:\n"
    "  path       array of positions, length N\n"
    "  histogram  array of bin counts\n"
    "  error      scatter across batches\n"
    "\n"
    "link_action(left, right, eps, m, hbar, V)\n"
    "  -> one term of eps * E\n"
    "\n"
    "direct_ground_density(...)\n"
    "  -> |phi|^2, value at 1, error\n"
    "\n"
    "metropolis_ground_state(...)\n"
    "  -> histogram, paths, correlator\n"
    "\n"
    "solve_schrodinger(V, m, hbar, grid)\n"
    "  -> wavefunction, density, energy\n"
    "\n"
    "INPUT   N, T, mass, hbar, V, x grid, samples\n"
    "OUTPUT  both densities, and |phi(1)|^2\n"
    "\n"
    "eps = T / N\n"
    "for each x:\n"
    "    average exp(-eps E) over bridges\n"
    "normalise\n"
    "compare with solve_schrodinger\n";

static const char* script_for(int page) {
    switch (page) {
    case PAGE_CONV: return kScriptConv;
    case PAGE_PATHS: return kScriptPaths;
    case PAGE_DIAG: return kScriptDiag;
    case PAGE_GAP: return kScriptGap;
    case PAGE_CLASS: return kScriptClass;
    case PAGE_HBAR: return kScriptHbar;
    case PAGE_POTS: return kScriptPots;
    case PAGE_DESIGN: return kScriptDesign;
    case PAGE_RESEARCH: return kScriptResearch;
    case PAGE_LIBRARY: return kScriptLibrary;
    default: return kScriptMile;
    }
}

static float g_scroll = 0;

static void draw_script(Rectangle r, const char* src) {
    DrawRectangleRounded(r, 0.03f, 6, Color{12, 16, 24, 255});
    DrawRectangleRoundedLines(r, 0.03f, 6, COL_LINE);
    DrawRectangle(r.x, r.y, 3, r.height, COL_TEAL);
    if (CheckCollisionPointRec(GetMousePosition(), r)) {
        g_scroll -= GetMouseWheelMove() * 28.0f;
        if (g_scroll < 0) g_scroll = 0;
    }
    BeginScissorMode((int)r.x, (int)r.y, (int)r.width, (int)r.height);
    float y = r.y + 12 - g_scroll;
    float limit = r.width - 28;
    const char* p = src;
    while (*p) {
        const char* nl = strchr(p, '\n');
        std::string raw(p, nl ? (size_t)(nl - p) : strlen(p));
        Color c = COL_INK;
        if (raw.rfind("INPUT", 0) == 0 || raw.rfind("OUTPUT", 0) == 0 || raw.rfind("DESIGN", 0) == 0) c = COL_GOLD;
        else if (raw.rfind("direct_", 0) == 0 || raw.rfind("metro", 0) == 0 || raw.rfind("solve_", 0) == 0 || raw.rfind("link_", 0) == 0 || raw.rfind("lattice_", 0) == 0 || raw.rfind("Potential", 0) == 0)
            c = COL_TEAL;
        if (raw.empty()) y += 18;
        size_t i = 0;
        while (i < raw.size()) {
            size_t take = raw.size() - i;
            while (take > 1 && text_width(g_code, raw.substr(i, take).c_str(), 14) > limit) take--;
            if (i + take < raw.size()) {
                size_t sp = raw.rfind(' ', i + take);
                if (sp != std::string::npos && sp > i + 6) take = sp - i;
            }
            std::string piece = raw.substr(i, take);
            if (y > r.y - 20 && y < r.y + r.height) text(g_code, piece.c_str(), r.x + 14, y, 14, c);
            y += 18;
            i += take;
            while (i < raw.size() && raw[i] == ' ') i++;
        }
        if (!nl) break;
        p = nl + 1;
    }
    EndScissorMode();
    float max_scroll = y + g_scroll - (r.y + r.height) + 20;
    if (max_scroll < 0) max_scroll = 0;
    if (g_scroll > max_scroll) g_scroll = max_scroll;
    if (max_scroll > 1.0f) {
        float track = r.height - 16;
        float thumb = std::max(24.0f, track * (r.height / (r.height + max_scroll)));
        float ty = r.y + 8 + (track - thumb) * (g_scroll / max_scroll);
        DrawRectangleRounded(Rectangle{r.x + r.width - 8, ty, 4, thumb}, 0.5f, 4, COL_TEAL);
    }
}

static void nudge(double& v, double lo, double hi, double step, int decimals) {
    if (decimals == 0) v = std::round(v);
    if (v < lo) v = lo;
    if (v > hi) v = hi;
}

static bool slider(Rectangle r, const char* name, double& v, double lo, double hi, double step, int decimals, bool logarithmic = false) {
    DrawRectangleRounded(r, 0.16f, 4, COL_INSET);
    text(g_ui, name, r.x + 8, r.y + 4, 13, COL_MUTED);
    char buf[32];
    if (decimals == 0) snprintf(buf, sizeof buf, "%.0f", v);
    else snprintf(buf, sizeof buf, "%.*f", decimals, v);
    float vw = text_width(g_ui, buf, 15);
    text(g_ui, buf, r.x + r.width - vw - 8, r.y + 4, 15, COL_INK);
    Rectangle track{r.x + 10, r.y + 32, r.width - 20, 4};
    DrawRectangleRounded(track, 0.5f, 4, COL_LINE);
    double llo = std::log(lo);
    double lhi = std::log(hi);
    float span = (float)(hi - lo);
    float t = 0;
    if (logarithmic && lhi > llo) t = (float)((std::log(std::max(v, lo)) - llo) / (lhi - llo));
    else if (span > 0) t = (float)((v - lo) / span);
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    float thumb = track.x + t * track.width;
    DrawRectangleRounded(Rectangle{track.x, track.y, thumb - track.x, 4}, 0.5f, 4, COL_TEAL);
    DrawCircleV(Vector2{thumb, track.y + 2}, 7, COL_INK);
    DrawCircleV(Vector2{thumb, track.y + 2}, 4, COL_TEAL);
    Rectangle hit{r.x, r.y + 18, r.width, r.height - 18};
    bool changed = false;
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), hit)) {
        float u = (GetMouseX() - track.x) / track.width;
        if (u < 0) u = 0;
        if (u > 1) u = 1;
        double nv = logarithmic ? std::exp(llo + (double)u * (lhi - llo)) : lo + (double)u * (hi - lo);
        nv = lo + step * std::round((nv - lo) / step);
        nudge(nv, lo, hi, step, decimals);
        if (std::fabs(nv - v) > 0.5 * step * 0.5) {
            v = nv;
            changed = true;
        }
    }
    return changed;
}

static void draw_controls(Rectangle area) {
    card(area);
    float sw = std::min(150.0f, (area.width - 188.0f) / 4.0f);
    if (sw < 78.0f) sw = 78.0f;
    float y = area.y + 6;
    float x0 = area.x + 10;
    bool moved = false;
    if (slider(Rectangle{x0, y, sw, 50}, "N", g_set_N, 4, 1000000, 1, 0, true)) moved = true;
    if (slider(Rectangle{x0 + (sw + 8), y, sw, 50}, "T", g_set_T, 1, 24, 0.5, 1)) moved = true;
    if (slider(Rectangle{x0 + 2 * (sw + 8), y, sw, 50}, "hbar", g_set_hbar, 0.2, 5, 0.1, 2)) moved = true;
    if (slider(Rectangle{x0 + 3 * (sw + 8), y, sw, 50}, "mass", g_set_mass, 0.25, 5, 0.25, 2)) moved = true;
    Rectangle play{area.x + area.width - 168, y + 8, 70, 32};
    if (button(play, g_play ? "Stop" : "Play", g_play)) g_play = !g_play;
    Rectangle run{area.x + area.width - 90, y + 8, 76, 32};
    if (button(run, g_busy ? "..." : "Run", !g_busy) && !g_busy) launch_run();
    if (moved) {
        g_play = false;
        refresh_exact();
        kick_live();
    }
    float y2 = y + 58;
    text(g_ui, "potential", area.x + 12, y2 + 6, 13, COL_MUTED);
    float x = area.x + 88;
    for (int i = 0; i < 5; ++i) {
        float bw = text_width(g_ui, kRunLabels[i], 14) + 22;
        if (x + bw > area.x + area.width - 8) break;
        if (button(Rectangle{x, y2, bw, 28}, kRunLabels[i], g_set_pot == i) && g_set_pot != i) {
            g_set_pot = i;
            g_play = false;
            refresh_exact();
            kick_live();
        }
        x += bw + 6;
    }
}

static void hist_bars(const Axis& a, const LargeHist& H) {
    if (!H.ok || H.x.size() < 2) return;
    double dx = H.x[1] - H.x[0];
    Rectangle clip = a.plot_rect();
    BeginScissorMode((int)clip.x, (int)clip.y, (int)clip.width, (int)clip.height);
    for (size_t i = 0; i < H.x.size(); ++i) {
        double xL = H.x[i] - 0.5 * dx;
        double xR = H.x[i] + 0.5 * dx;
        if (xR < a.x0 || xL > a.x1) continue;
        Vector2 lo = a.map(xL, 0.0);
        Vector2 hi = a.map(xR, H.psi2[i]);
        float top = std::min(lo.y, hi.y);
        float bot = std::max(lo.y, hi.y);
        DrawRectangleRec(Rectangle{lo.x + 1, top, std::max(1.0f, hi.x - lo.x - 1), std::max(0.0f, bot - top)}, Color{32, 148, 140, 150});
    }
    EndScissorMode();
}

static void draw_milestone(Rectangle area, const Pack& p) {
    card(area);
    bool ho = g_live.ready ? (g_live.pot == 0) : (p.run_pot && strcmp(p.run_pot, "ho") == 0);
    int N = g_live.ready ? g_live.N : p.N;
    double T = g_live.ready ? g_live.T : p.T;
    double mass = g_live.ready ? g_live.mass : p.mass;
    double hbar = g_live.ready ? g_live.hbar : p.hbar_run;
    const Spectrum& se = g_live.ready ? g_live.se : p.schrodinger;
    const HoLattice& lat = g_live.ready ? g_live.lattice : p.lattice;
    bool matched = sample_matches(p);
    bool large = N > 48;
    bool hmatch = hist_matches(p);
    text(g_ui, large ? "Histogram of the slices" : (ho ? "Ground state from x = 0 to 2" : "Ground-state density"), area.x + 16, area.y + 12, 22, COL_INK);
    char sub[240];
    const char* tail = "";
    if (large && hmatch) tail = "    data/histogram.csv";
    else if (large) tail = "    Run writes data/histogram.csv";
    else if (!matched) tail = "    sampling";
    snprintf(sub, sizeof sub, "N = %d    T = %g    eps = %.4g    hbar = %g    mass = %g    %s%s", N, T, T / std::max(N, 1), hbar, mass, kRunPots[g_live.ready ? g_live.pot : 0], tail);
    text(g_ui, sub, area.x + 16, area.y + 40, 15, (matched || hmatch) ? COL_MUTED : COL_GOLD);
    double x0 = 0, x1 = 2;
    if (large) {
        if (ho) {
            x0 = -6;
            x1 = 6;
        } else {
            const char* pname = kRunPots[g_live.ready ? g_live.pot : 0];
            pot_window(pname, x0, x1);
        }
    } else if (!ho) {
        if (matched && !p.direct.x.empty()) {
            x0 = p.direct.x.front();
            x1 = p.direct.x.back();
        } else if (!se.x.empty()) {
            x0 = se.x.front();
            x1 = se.x.back();
        }
    }
    double ymax = 0.05;
    auto bump = [&](const std::vector<double>& xs, const std::vector<double>& ys) {
        size_t n = std::min(xs.size(), ys.size());
        for (size_t i = 0; i < n; ++i) {
            if (xs[i] < x0 || xs[i] > x1) continue;
            if (ys[i] > ymax) ymax = ys[i];
        }
    };
    if (matched) bump(p.direct.x, p.direct.psi2);
    if (hmatch) bump(p.hist.x, p.hist.psi2);
    bump(se.x, se.psi2);
    if (ho && lat.ok) {
        double peak = gaussian_psi2(0.0, lat.alpha);
        if (peak > ymax) ymax = peak;
    }
    Axis a = make_axis(Rectangle{area.x + 8, area.y + 64, area.width - 16, area.height - 80}, x0, x1, 0, ymax * 1.18);
    draw_frame(a, "x", "|phi(x)|^2");
    if (hmatch) hist_bars(a, p.hist);
    if (matched) band(a, p.direct.x, p.direct.psi2, p.direct.err, COL_TEAL_DIM);
    std::vector<double> lx, ly, sx, sy, ax, ay;
    if (ho) {
        for (int i = 0; i <= 200; ++i) {
            double x = x0 + (x1 - x0) * i / 200.0;
            ax.push_back(x);
            ay.push_back(ho_ground_psi2(x, mass, hbar, 1.0));
            if (lat.ok) {
                lx.push_back(x);
                ly.push_back(gaussian_psi2(x, lat.alpha));
            }
        }
    }
    for (size_t i = 0; i < se.x.size(); ++i) {
        sx.push_back(se.x[i]);
        sy.push_back(se.psi2[i]);
    }
    stroke(a, lx, ly, COL_NAVY, 2.0f, false);
    stroke(a, sx, sy, COL_INK, 2.2f, false);
    stroke(a, ax, ay, COL_COPPER, 1.6f, true);
    if (matched) stroke(a, p.direct.x, p.direct.psi2, COL_TEAL, 1.6f, false);
    if (x0 <= 1.0 && x1 >= 1.0) {
        double ydot = matched ? p.direct.psi2_at_1 : (lat.ok ? gaussian_psi2(1.0, lat.alpha) : (se.ok ? interp(se.x, se.psi2, 1.0) : 0.0));
        Vector2 dot = a.map(1.0, ydot);
        DrawCircleV(dot, 5, COL_TEAL);
        DrawCircleLines((int)dot.x, (int)dot.y, 5, COL_INK);
    }
    if (large && ho) {
        const char* names[] = {"histogram", "Schrodinger", "analytic"};
        Color cols[] = {COL_TEAL, COL_INK, COL_COPPER};
        legend(a.plot_rect(), names, cols, 3);
    } else if (large) {
        const char* names[] = {"histogram", "Schrodinger"};
        Color cols[] = {COL_TEAL, COL_INK};
        legend(a.plot_rect(), names, cols, 2);
    } else if (ho) {
        const char* names[] = {"Monte Carlo band", "exact lattice", "Schrodinger", "analytic", "direct paths"};
        Color cols[] = {COL_TEAL_DIM, COL_NAVY, COL_INK, COL_COPPER, COL_TEAL};
        legend(a.plot_rect(), names, cols, 5);
    } else {
        const char* names[] = {"Monte Carlo band", "Schrodinger", "direct paths"};
        Color cols[] = {COL_TEAL_DIM, COL_INK, COL_TEAL};
        legend(a.plot_rect(), names, cols, 3);
    }
}

static void draw_convergence(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Percent above the true ground state at x = 1", area.x + 16, area.y + 12, 20, COL_INK);
    float grid_x = area.x + 78;
    float grid_y = area.y + 86;
    float cw = (area.width - 110) / 6.0f;
    float ch = (area.height - 160) / 7.0f;
    if (ch > 78) ch = 78;
    for (int j = 0; j < 6; ++j) {
        char b[16];
        snprintf(b, sizeof b, "%g", p.Ts[j]);
        text(g_ui, b, grid_x + j * cw + cw * 0.3f, grid_y - 22, 14, COL_MUTED);
    }
    text(g_ui, "T", grid_x + 3 * cw, area.y + 46, 14, COL_MUTED);
    double ground = p.analytic1 > 0 ? p.analytic1 : 0.2075537487;
    for (int i = 0; i < 7; ++i) {
        char b[16];
        snprintf(b, sizeof b, "%d", p.Ns[6 - i]);
        text(g_ui, b, area.x + 28, grid_y + i * ch + ch * 0.3f, 14, COL_MUTED);
        for (int j = 0; j < 6; ++j) {
            int row = 6 - i;
            double rel = 100.0 * (p.psi1[row][j] - ground) / ground;
            float shade = (float)std::min(1.0, std::max(0.0, rel / 12.0));
            Color cell{
                (unsigned char)(40 + shade * 170),
                (unsigned char)(36 + (1.0f - shade) * 40),
                (unsigned char)(32 + (1.0f - shade) * 20),
                255
            };
            Rectangle rc{grid_x + j * cw + 2, grid_y + i * ch + 2, cw - 4, ch - 4};
            DrawRectangleRounded(rc, 0.15f, 4, cell);
            snprintf(b, sizeof b, "%.1f", rel);
            Color tc = rel > 6 ? COL_INK : COL_MUTED;
            text(g_ui, b, rc.x + 8, rc.y + rc.height * 0.3f, 15, tc);
            if (p.Ns[row] == 7 && std::fabs(p.Ts[j] - 4) < 0.01) DrawRectangleLinesEx(rc, 2, COL_TEAL);
        }
    }
    text(g_ui, "N", area.x + 24, area.y + 58, 14, COL_MUTED);
    text(g_ui, "Teal box is the milestone, N = 7, T = 4.", area.x + 16, area.y + area.height - 28, 15, COL_MUTED);
}

static Texture2D g_world{};
static int g_world_gen = -1;

static void ensure_world(const Pack& p) {
    if (g_world_gen == p.generation && g_world.id != 0) return;
    if (p.metro.world.empty() || p.metro.world_bins < 2 || p.metro.N < 2) return;
    if (g_world.id) UnloadTexture(g_world);
    int bins = p.metro.world_bins;
    int N = p.metro.N;
    Image img = GenImageColor(N, bins, COL_BG);
    double maxd = 1e-12;
    for (double v : p.metro.world) if (v > maxd) maxd = v;
    for (int j = 0; j < N; ++j) {
        for (int b = 0; b < bins; ++b) {
            double d = p.metro.world[(size_t)j * bins + b] / maxd;
            unsigned char r = (unsigned char)(16 + d * 220);
            unsigned char g = (unsigned char)(24 + d * 150);
            unsigned char bl = (unsigned char)(32 + (1.0 - d) * 40 + d * 70);
            ImageDrawPixel(&img, j, bins - 1 - b, Color{r, g, bl, 255});
        }
    }
    g_world = LoadTextureFromImage(img);
    UnloadImage(img);
    SetTextureFilter(g_world, TEXTURE_FILTER_BILINEAR);
    g_world_gen = p.generation;
}

static void draw_paths(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Imaginary-time paths", area.x + 16, area.y + 12, 22, COL_INK);
    ensure_world(p);
    Rectangle heat{area.x + 16, area.y + 48, area.width - 32, area.height * 0.52f};
    DrawRectangleRec(heat, COL_INSET);
    if (g_world.id) {
        DrawTexturePro(g_world, Rectangle{0, 0, (float)g_world.width, (float)g_world.height}, heat, Vector2{0, 0}, 0, WHITE);
    }
    float y1 = heat.y + heat.height * (float)((3.0 - p.metro.world_xmin) / (p.metro.world_xmax - p.metro.world_xmin));
    float y2 = heat.y + heat.height * (float)((-3.0 - p.metro.world_xmin) / (p.metro.world_xmax - p.metro.world_xmin));
    (void)y2;
    auto y_of = [&](double x) {
        double t = (x - p.metro.world_xmin) / (p.metro.world_xmax - p.metro.world_xmin);
        return heat.y + heat.height * (float)(1.0 - t);
    };
    DrawLineEx(Vector2{heat.x, y_of(1)}, Vector2{heat.x + heat.width, y_of(1)}, 1, Color{204, 92, 58, 180});
    DrawLineEx(Vector2{heat.x, y_of(-1)}, Vector2{heat.x + heat.width, y_of(-1)}, 1, Color{204, 92, 58, 180});
    for (size_t i = 0; i < p.metro.paths.size() && i < 8; ++i) {
        const auto& path = p.metro.paths[i];
        for (size_t j = 1; j < path.size(); ++j) {
            float x0 = heat.x + heat.width * (float)(j - 1) / (float)path.size();
            float x1 = heat.x + heat.width * (float)j / (float)path.size();
            DrawLineEx(Vector2{x0, y_of(path[j - 1])}, Vector2{x1, y_of(path[j])}, 1.1f, Color{246, 231, 193, 150});
        }
    }
    text(g_ui, "tau", heat.x + heat.width - 36, heat.y + heat.height - 22, 14, COL_MUTED);
    text(g_ui, "copper lines: classical turning points", heat.x + 8, heat.y + 8, 14, COL_INK);

    Rectangle br{area.x + 16, heat.y + heat.height + 16, area.width - 32, area.y + area.height - (heat.y + heat.height + 28)};
    DrawRectangleRec(br, COL_INSET);
    text(g_ui, "Milestone bridges pinned at x = 1, coloured by weight", br.x + 8, br.y + 6, 14, COL_MUTED);
    if (p.bridges.size() > 1 && !p.bridges[0].x.empty()) {
        double lo = 1e300, hi = -1e300;
        for (const auto& b : p.bridges) {
            if (!(b.weight > 0)) continue;
            double lw = std::log(b.weight);
            lo = std::min(lo, lw);
            hi = std::max(hi, lw);
        }
        if (!(hi > lo)) hi = lo + 1;
        int n = (int)p.bridges[0].x.size();
        for (const auto& b : p.bridges) {
            float shade = (b.weight > 0) ? (float)((std::log(b.weight) - lo) / (hi - lo)) : 0;
            Color c{(unsigned char)(70 + shade * 180), (unsigned char)(80 + shade * 120), (unsigned char)(90 - shade * 40), 230};
            for (int j = 1; j < n; ++j) {
                float x0 = br.x + 10 + (br.width - 20) * (float)(j - 1) / (float)(n - 1);
                float x1 = br.x + 10 + (br.width - 20) * (float)j / (float)(n - 1);
                float y0 = br.y + br.height - 16 - (br.height - 40) * (float)((b.x[j - 1] + 1.5) / 5.0);
                float y1 = br.y + br.height - 16 - (br.height - 40) * (float)((b.x[j] + 1.5) / 5.0);
                DrawLineEx(Vector2{x0, y0}, Vector2{x1, y1}, 1.3f, c);
            }
        }
    }
}

static void draw_diagnostics(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Chain, correlator, autocorrelation", area.x + 16, area.y + 12, 22, COL_INK);
    float w = (area.width - 36) * 0.5f;
    float h = (area.height - 88) * 0.5f;
    Rectangle boxes[4] = {
        {area.x + 12, area.y + 48, w, h},
        {area.x + 24 + w, area.y + 48, w, h},
        {area.x + 12, area.y + 60 + h, w, h},
        {area.x + 24 + w, area.y + 60 + h, w, h}
    };
    const MetroResult& m = p.metro;
    if (m.trace_x2.size() > 4) {
        double xmin = m.sweep.empty() ? 0 : m.sweep.front();
        double xmax = m.sweep.empty() ? (double)m.trace_x2.size() : m.sweep.back();
        Axis a = make_axis(boxes[0], xmin, xmax, 0.2, 0.9);
        draw_frame(a, "sweep", "<x^2>");
        stroke(a, m.sweep, m.trace_x2, COL_TEAL, 1.4f, false);
        std::vector<double> hx{xmin, xmax};
        std::vector<double> hy{0.5, 0.5};
        stroke(a, hx, hy, COL_COPPER, 1.2f, true);
    }
    if (m.trace_acc.size() > 4) {
        double xmin = m.sweep.empty() ? 0 : m.sweep.front();
        double xmax = m.sweep.empty() ? (double)m.trace_acc.size() : m.sweep.back();
        Axis a = make_axis(boxes[1], xmin, xmax, 0, 1);
        draw_frame(a, "sweep", "acceptance");
        stroke(a, m.sweep, m.trace_acc, COL_GOLD, 1.4f, false);
    }
    if (m.corr_tau.size() > 2) {
        double tmax = m.corr_tau.back();
        double ymax = 0.2;
        for (double v : m.corr_xx) ymax = std::max(ymax, v);
        Axis a = make_axis(boxes[2], 0, tmax, 0, ymax * 1.15);
        draw_frame(a, "tau", "<x(tau) x(0)>");
        std::vector<double> tx, ty;
        for (int i = 0; i <= 80; ++i) {
            double tau = tmax * i / 80.0;
            double s = std::sinh(0.5 * m.T);
            tx.push_back(tau);
            ty.push_back((s > 0) ? 0.5 * std::cosh(0.5 * m.T - tau) / s : 0);
        }
        band(a, m.corr_tau, m.corr_xx, m.corr_err, COL_TEAL_DIM);
        stroke(a, m.corr_tau, m.corr_xx, COL_TEAL, 2.0f, false);
        stroke(a, tx, ty, COL_INK, 1.4f, true);
    }
    if (p.acf.size() > 2) {
        std::vector<double> lx, ly;
        double ymax = 0.05;
        for (size_t i = 0; i < p.acf.size(); ++i) {
            lx.push_back((double)i);
            ly.push_back(p.acf[i]);
            ymax = std::max(ymax, std::fabs(p.acf[i]));
        }
        Axis a = make_axis(boxes[3], 0, (double)p.acf.size(), -0.2, std::max(1.05, ymax));
        draw_frame(a, "lag (sweeps)", "autocorrelation");
        stroke(a, lx, ly, COL_COPPER, 2.0f, false);
    }
    char buf[96];
    snprintf(buf, sizeof buf, "gap fit %.2f    true gap 1", p.gap);
    text(g_ui, buf, area.x + 16, area.y + area.height - 24, 15, COL_MUTED);
}

static void draw_classical(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Classical particle against the quantum one", area.x + 16, area.y + 12, 22, COL_INK);
    Axis a = make_axis(Rectangle{area.x + 8, area.y + 48, area.width - 16, area.height - 70}, -2.4, 2.4, 0, 1.15);
    draw_frame(a, "x", "probability density");
    BeginScissorMode((int)a.plot_rect().x, (int)a.plot_rect().y, (int)a.plot_rect().width, (int)a.plot_rect().height);
    Vector2 left = a.map(-2.4, 0);
    Vector2 left2 = a.map(-1, 1.15);
    DrawRectangleRec(Rectangle{left.x, a.plot_rect().y, left2.x - left.x, a.plot_rect().height}, Color{80, 62, 48, 90});
    Vector2 right = a.map(1, 0);
    DrawRectangleRec(Rectangle{right.x, a.plot_rect().y, a.map(2.4, 0).x - right.x, a.plot_rect().height}, Color{80, 62, 48, 90});
    EndScissorMode();
    std::vector<double> cx, cy, qx, qy, ax, ay;
    for (int i = 0; i < 400; ++i) {
        double x = -0.985 + 1.97 * i / 399.0;
        cx.push_back(x);
        cy.push_back(1.0 / (3.141592653589793 * std::sqrt(1.0 - x * x)));
        double z = -2.4 + 4.8 * i / 399.0;
        ax.push_back(z);
        ay.push_back(std::exp(-z * z) / std::sqrt(3.141592653589793));
    }
    if (p.metro.ok) band(a, p.metro.x, p.metro.psi2, p.metro.err, COL_TEAL_DIM);
    stroke(a, cx, cy, COL_COPPER, 2.0f, false);
    if (p.metro.ok) stroke(a, p.metro.x, p.metro.psi2, COL_TEAL, 2.0f, false);
    stroke(a, ax, ay, COL_INK, 1.5f, true);
    const char* names[] = {"classical", "path integral", "analytic"};
    Color cols[] = {COL_COPPER, COL_TEAL, COL_INK};
    legend(a.plot_rect(), names, cols, 3);
    char buf[160];
    snprintf(buf, sizeof buf, "Chance past |x| = 1:   path integral %.1f%%     exact ground state %.1f%%", 100.0 * p.forbid, 100.0 * std::erfc(1.0));
    text(g_ui, buf, area.x + 16, area.y + area.height - 26, 15, COL_MUTED);
}

static void draw_hbar(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Turn hbar down and the packet gets narrower", area.x + 16, area.y + 12, 22, COL_INK);
    Axis a = make_axis(Rectangle{area.x + 8, area.y + 48, area.width - 16, area.height - 64}, -6.5, 6.5, -0.15, 6.3);
    draw_frame(a, "x", "density, stacked");
    Color ridge[6] = {
        Color{120, 42, 36, 230}, Color{196, 92, 58, 230}, Color{180, 140, 60, 230},
        Color{32, 130, 124, 230}, Color{70, 100, 150, 230}, Color{110, 70, 120, 230}
    };
    for (int i = 0; i < 6; ++i) {
        if (!p.hbar[i].ok) continue;
        double base = (5 - i) * 1.05;
        std::vector<double> ys = p.hbar[i].psi2;
        std::vector<double> an;
        an.reserve(p.hbar[i].x.size());
        for (double x : p.hbar[i].x) an.push_back(ho_ground_psi2(x, 1.0, p.hbars[i], 1.0));
        for (double& y : ys) y += base;
        for (double& y : an) y += base;
        band(a, p.hbar[i].x, ys, p.hbar[i].err, Color{ridge[i].r, ridge[i].g, ridge[i].b, 50});
        stroke(a, p.hbar[i].x, ys, ridge[i], 2.0f, false);
        stroke(a, p.hbar[i].x, an, COL_INK, 1.0f, false);
        char lab[32];
        snprintf(lab, sizeof lab, "hbar %g", p.hbars[i]);
        Vector2 pos = a.map(-6.3, base + 0.15);
        text(g_ui, lab, pos.x, pos.y - 16, 15, ridge[i]);
    }
}

static int g_pot = 1;

static void draw_pots(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Same engine, other potentials", area.x + 16, area.y + 12, 22, COL_INK);
    const char* labels[4] = {"Anharmonic", "Double well", "Morse", "Quartic"};
    for (int i = 0; i < 4; ++i) {
        Rectangle b{area.x + 16 + i * 150.0f, area.y + 46, 140, 30};
        if (button(b, labels[i], g_pot == i)) g_pot = i;
    }
    int i = g_pot;
    if (i < 0 || i > 3 || !p.pot_metro[i].ok) return;
    const MetroResult& m = p.pot_metro[i];
    const Spectrum& se = p.pot_se[i];
    Axis a = make_axis(Rectangle{area.x + 8, area.y + 88, area.width - 16, area.height - 120}, m.world_xmin, m.world_xmax, 0, 0);
    double ymax = 0.05;
    for (double v : m.psi2) ymax = std::max(ymax, v);
    for (double v : se.psi2) ymax = std::max(ymax, v);
    a.y1 = ymax * 1.25;
    a.x0 = m.x.empty() ? -3 : m.x.front();
    a.x1 = m.x.empty() ? 3 : m.x.back();
    draw_frame(a, "x", "|phi(x)|^2");
    Potential pot = parse_potential(p.pot_name[i]);
    std::vector<double> vx, vy;
    double vmax = 1e-6;
    for (double x : m.x) vmax = std::max(vmax, pot(x));
    for (double x : m.x) {
        vx.push_back(x);
        vy.push_back(pot(x) / vmax * a.y1 * 0.92);
    }
    stroke(a, vx, vy, Color{196, 160, 120, 180}, 1.4f, false);
    if (se.ok) {
        double e = se.energy / vmax * a.y1 * 0.92;
        Vector2 e0 = a.map(a.x0, e);
        Vector2 e1 = a.map(a.x1, e);
        DrawLineEx(e0, e1, 1.2f, COL_COPPER);
        stroke(a, se.x, se.psi2, COL_INK, 2.0f, false);
    }
    band(a, m.x, m.psi2, m.err, COL_TEAL_DIM);
    stroke(a, m.x, m.psi2, COL_TEAL, 2.0f, false);
    const char* names[] = {"potential, scaled", "energy", "Schrodinger", "path integral"};
    Color cols[] = {Color{196, 160, 120, 255}, COL_COPPER, COL_INK, COL_TEAL};
    legend(a.plot_rect(), names, cols, 4);
    char buf[128];
    if (se.ok) snprintf(buf, sizeof buf, "Ground energy  %.4f     add another case in potential.hpp", se.energy);
    else snprintf(buf, sizeof buf, "Schrodinger solve failed for this potential");
    text(g_ui, buf, area.x + 16, area.y + area.height - 26, 15, COL_MUTED);
}

static void draw_gap(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "The splitting, and why the chain misses it", area.x + 18, area.y + 14, 22, COL_INK);
    text(g_ui, "Same Hamiltonian. Spectrum on the left of each pair, paths on the right.", area.x + 18, area.y + 42, 15, COL_MUTED);
    float mid = area.x + area.width * 0.5f;
    auto stat = [&](float x, float y, const char* k, const char* v, Color c) {
        text(g_ui, k, x, y, 14, COL_MUTED);
        text(g_ui, v, x + 118, y, 15, c);
    };
    char b0[64], b1[64], b2[64], b3[64];
    text(g_ui, "Oscillator", area.x + 18, area.y + 74, 16, COL_TEAL);
    snprintf(b0, sizeof b0, "%.3f", p.ho_e0);
    snprintf(b1, sizeof b1, "%.3f", p.ho_e1);
    snprintf(b2, sizeof b2, "%.3f", p.ho_e1 - p.ho_e0);
    snprintf(b3, sizeof b3, "%.2f  per path", p.ho_cross);
    stat(area.x + 18, area.y + 100, "ground E0", b0, COL_INK);
    stat(area.x + 18, area.y + 122, "excited E1", b1, COL_INK);
    stat(area.x + 18, area.y + 144, "gap E1-E0", b2, COL_GOLD);
    stat(area.x + 18, area.y + 166, "zero crossings", b3, COL_INK);
    char energy[64];
    snprintf(energy, sizeof energy, "%.2f ± %.2f", p.ho_energy, p.ho_energy_err);
    stat(area.x + 18, area.y + 188, "path energy", energy, COL_INK);
    text(g_ui, "against 0.50", area.x + 136 + text_width(g_ui, energy, 15) + 8, area.y + 188, 14, COL_MUTED);
    text(g_ui, "Double well", mid + 8, area.y + 74, 16, COL_COPPER);
    snprintf(b0, sizeof b0, "%.3f", p.dw_e0);
    snprintf(b1, sizeof b1, "%.3f", p.dw_e1);
    snprintf(b2, sizeof b2, "%.3f", p.dw_e1 - p.dw_e0);
    snprintf(b3, sizeof b3, "%.2f  per path", p.dw_cross);
    stat(mid + 8, area.y + 100, "ground E0", b0, COL_INK);
    stat(mid + 8, area.y + 122, "excited E1", b1, COL_INK);
    stat(mid + 8, area.y + 144, "gap E1-E0", b2, COL_GOLD);
    stat(mid + 8, area.y + 166, "zero crossings", b3, COL_INK);

    float plot_top = area.y + 236;
    Rectangle plot{area.x + 58, plot_top, area.width - 86, area.height - (plot_top - area.y) - 36};
    if (plot.height < 40) return;
    DrawRectangleRec(plot, COL_INSET);
    text(g_ui, "double-well gap against hbar", plot.x, plot.y - 20, 14, COL_MUTED);
    double ymin = -5.0;
    double ymax = 0.4;
    auto X = [&](double h) { return plot.x + (float)((h - 0.40) / (2.30 - 0.40)) * plot.width; };
    auto Y = [&](double g) {
        double ly = std::log10(std::max(g, 1e-6));
        float t = (float)((ly - ymin) / (ymax - ymin));
        return plot.y + plot.height - t * plot.height;
    };
    for (int i = 0; i < 5; ++i) {
        float y = plot.y + plot.height * (float)i / 4.0f;
        DrawLineEx(Vector2{plot.x, y}, Vector2{plot.x + plot.width, y}, 1, COL_LINE);
    }
    for (int i = 0; i < 6; ++i) {
        Vector2 a{X(p.split_h[i]), Y(p.split_gap[i])};
        Vector2 b{X(p.split_h[i + 1]), Y(p.split_gap[i + 1])};
        DrawLineEx(a, b, 2.4f, COL_COPPER);
        DrawCircleV(a, 4.0f, COL_COPPER);
    }
    DrawCircleV(Vector2{X(p.split_h[6]), Y(p.split_gap[6])}, 4.0f, COL_COPPER);
    float hx = X(1.0);
    DrawLineEx(Vector2{hx, plot.y}, Vector2{hx, plot.y + plot.height}, 1.2f, COL_TEAL);
    text(g_ui, "hbar = 1", hx + 6, plot.y + 8, 13, COL_TEAL);
    const char* ylabs[] = {"1", "0.1", "0.01", "10^-3", "10^-4"};
    for (int i = 0; i < 5; ++i) {
        float y = Y(std::pow(10.0, -i));
        text(g_ui, ylabs[i], area.x + 12, y - 7, 13, COL_MUTED);
    }
    text(g_ui, "0.5", X(0.5) - 8, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "1", X(1.0) - 4, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "2", X(2.0) - 4, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "hbar", plot.x + plot.width - 36, plot.y + plot.height + 6, 13, COL_MUTED);
}

static void draw_research(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "Two checks the milestone does not ask for", area.x + 18, area.y + 14, 22, COL_INK);
    text(g_ui, "Classical tail, and the gap where there is no formula.", area.x + 18, area.y + 42, 15, COL_MUTED);
    char row[160];
    snprintf(row, sizeof row, "gap    oscillator 1.000    anharmonic %.3f    Morse %.3f    quartic %.3f", p.an_gap, p.mo_gap, p.qu_gap);
    text(g_ui, row, area.x + 18, area.y + 68, 15, COL_INK);
    int i1 = 0;
    for (int i = 1; i < 7; ++i) if (std::fabs(p.split_h[i] - 1.0) < std::fabs(p.split_h[i1] - 1.0)) i1 = i;
    snprintf(row, sizeof row, "double well    hbar 1   barrier %.0f%%   gap %.3f       hbar 0.45   barrier %.0f%%   gap %.5f", 100.0 * p.res_bar[i1], p.split_gap[i1], 100.0 * p.res_bar[0], p.split_gap[0]);
    text(g_ui, row, area.x + 18, area.y + 90, 15, COL_MUTED);

    float plot_top = area.y + 156;
    Rectangle plot{area.x + 64, plot_top, area.width - 96, area.height - (plot_top - area.y) - 72};
    if (plot.height < 40) return;
    DrawRectangleRec(plot, COL_INSET);
    text(g_ui, "chance of |x| > 1", plot.x, plot.y - 20, 14, COL_MUTED);
    double ymin = 0.0, ymax = 0.55;
    auto X = [&](double h) { return plot.x + (float)((h - 0.30) / (3.80 - 0.30)) * plot.width; };
    auto Y = [&](double v) {
        float t = (float)((v - ymin) / (ymax - ymin));
        return plot.y + plot.height - t * plot.height;
    };
    for (int i = 0; i < 4; ++i) {
        float y = plot.y + plot.height * (float)i / 4.0f;
        DrawLineEx(Vector2{plot.x, y}, Vector2{plot.x + plot.width, y}, 1, COL_LINE);
    }
    float yturn = Y(ho_outside_turning(1.0));
    DrawLineEx(Vector2{plot.x, yturn}, Vector2{plot.x + plot.width, yturn}, 1.4f, COL_NAVY);
    text(g_ui, "turning point  15.7%", plot.x + plot.width - 150, yturn + 6, 13, COL_NAVY);
    for (int i = 0; i < 7; ++i) {
        DrawLineEx(Vector2{X(p.res_h[i]), Y(p.res_exact[i])}, Vector2{X(p.res_h[i + 1]), Y(p.res_exact[i + 1])}, 2.4f, COL_COPPER);
    }
    for (int i = 0; i < 6; ++i) DrawCircleV(Vector2{X(p.hbars[i]), Y(p.res_mc[i])}, 4.5f, COL_TEAL);
    text(g_ui, "0", X(0.35) - 4, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "1", X(1.0) - 4, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "2", X(2.0) - 4, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "3", X(3.0) - 4, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "hbar", plot.x + plot.width - 32, plot.y + plot.height + 6, 13, COL_MUTED);
    text(g_ui, "Copper is exact. Teal dots are the Metropolis scan.", area.x + 18, area.y + area.height - 46, 14, COL_MUTED);
    snprintf(row, sizeof row, "A short solve %.1f ms, a short chain %.1f ms, the direct sum %.1f ms.", p.ms_se, p.ms_metro, p.ms_direct);
    text(g_ui, row, area.x + 18, area.y + area.height - 26, 14, COL_MUTED);
}

static void draw_library(Rectangle area, const Pack& p) {
    (void)p;
    card(area);
    text(g_ui, "The library, and what to build next", area.x + 18, area.y + 14, 22, COL_INK);
    text(g_ui, "The window is one program. The functions live in cpp/feynman and can be included on their own.", area.x + 18, area.y + 44, 15, COL_MUTED);
    float colw = (area.width - 48) / 2.0f;
    Rectangle left{area.x + 16, area.y + 78, colw, area.height - 110};
    Rectangle right{area.x + 32 + colw, area.y + 78, colw, area.height - 110};
    DrawRectangleRounded(left, 0.04f, 6, COL_INSET);
    DrawRectangleRounded(right, 0.04f, 6, COL_INSET);
    text(g_ui, "Documents", left.x + 14, left.y + 12, 16, COL_GOLD);
    const char* docs[] = {
        "WRITEUP.md", "Milestone, design, the figures.",
        "THEORY.md", "Why the chain misses the splitting.",
        "README.md", "Build, suite, and this window.",
        "feynman.hpp", "One include for a new simulation."
    };
    float y = left.y + 42;
    for (int i = 0; i < 4; ++i) {
        text(g_ui, docs[i * 2], left.x + 14, y, 16, COL_INK);
        text(g_ui, docs[i * 2 + 1], left.x + 14, y + 20, 14, COL_MUTED);
        y += 48;
    }
    text(g_ui, "So the next run is cheaper", right.x + 14, right.y + 12, 16, COL_GOLD);
    const char* steps[] = {
        "Include the header. Do not copy the window.",
        "Histogram one chain. Do not bridge every x.",
        "Keep the gallery. A new question is a new file.",
        "On the double well, add a kink move.",
        "Split chains across threads, as the direct sum already splits x."
    };
    y = right.y + 44;
    for (int i = 0; i < 5; ++i) {
        char n[8];
        snprintf(n, sizeof n, "%d", i + 1);
        text(g_ui, n, right.x + 14, y, 16, COL_TEAL);
        text(g_ui, steps[i], right.x + 36, y, 15, COL_INK);
        y += 36;
    }
    text(g_ui, "The long chains in the gallery are the expensive part. Extra local sweeps on the double well do not buy the gap.", area.x + 18, area.y + area.height - 26, 14, COL_MUTED);
}

static void draw_design(Rectangle area, const Pack& p) {
    card(area);
    text(g_ui, "What the brief asks, and where it is", area.x + 16, area.y + 12, 22, COL_INK);
    struct Row { const char* item; const char* where; bool ok; };
    char mile[80], metro[80];
    snprintf(mile, sizeof mile, "direct |phi(1)|^2 = %.5f +/- %.5f", p.direct.psi2_at_1, p.direct.err_at_1);
    snprintf(metro, sizeof metro, "long chain |phi(1)|^2 = %.4f,  <x^2> = %.3f", p.metro.psi2_at_1, p.metro.mean_x2);
    Row rows[5] = {
        {"Design: functions, inputs, outputs, data structures", "this page and the script", true},
        {"Plot |phi|^2 from x = 0 to 2, both methods", "Milestone", p.direct.ok && p.schrodinger.ok},
        {mile, "Milestone, the dot at x = 1", p.direct.ok},
        {"Metropolis reproduces the ground state", metro, p.metro.ok},
        {"Room to add a potential or a page", "Potentials, kPotOrder, kPages", true}
    };
    float y = area.y + 56;
    for (int i = 0; i < 5; ++i) {
        DrawRectangleRounded(Rectangle{area.x + 16, y, area.width - 32, 62}, 0.12f, 6, COL_INSET);
        DrawRectangle(area.x + 16, y, 5, 62, rows[i].ok ? COL_TEAL : COL_COPPER);
        text(g_ui, rows[i].ok ? "done" : "open", area.x + 32, y + 8, 14, rows[i].ok ? COL_TEAL : COL_COPPER);
        text(g_ui, rows[i].item, area.x + 84, y + 8, 16, COL_INK);
        text(g_ui, rows[i].where, area.x + 84, y + 32, 14, COL_MUTED);
        y += 72;
    }
    text(g_ui, "Schrodinger at x = 1 matches the analytic 0.20755. The path integral on N = 7, T = 4 is a bit higher because T is short.", area.x + 16, area.y + area.height - 28, 14, COL_MUTED);
}

static void draw_numbers(Rectangle r, int page, const Pack& p) {
    card(r);
    text(g_ui, "At x = 1", r.x + 14, r.y + 12, 16, COL_MUTED);
    char big[64];
    bool live_num = page == PAGE_MILE && g_live.ready;
    if (live_num && hist_matches(p)) snprintf(big, sizeof big, "%.5f", interp(p.hist.x, p.hist.psi2, 1.0));
    else if (live_num && sample_matches(p)) snprintf(big, sizeof big, "%.5f", p.direct.psi2_at_1);
    else if (live_num && g_live.lattice.ok) snprintf(big, sizeof big, "%.5f", gaussian_psi2(1.0, g_live.lattice.alpha));
    else if (live_num && g_live.se.ok) snprintf(big, sizeof big, "%.5f", interp(g_live.se.x, g_live.se.psi2, 1.0));
    else if (page == PAGE_MILE || page == PAGE_DESIGN) snprintf(big, sizeof big, "%.5f", p.direct.psi2_at_1);
    else if (page == PAGE_PATHS) snprintf(big, sizeof big, "%.4f", p.metro.psi2_at_1);
    else if (page == PAGE_CLASS) snprintf(big, sizeof big, "%.1f%%", 100.0 * p.forbid);
    else if (page == PAGE_POTS && p.pot_se[g_pot].ok) snprintf(big, sizeof big, "%.4f", p.pot_se[g_pot].energy);
    else snprintf(big, sizeof big, "%.5f", p.analytic1);
    text(g_ui, big, r.x + 14, r.y + 34, 32, COL_INK);
    const char* cap = "path integral";
    if (live_num && hist_matches(p)) cap = "histogram at x = 1";
    if (live_num && !sample_matches(p) && !hist_matches(p) && g_live.lattice.ok) cap = "exact lattice, this N";
    if (live_num && !sample_matches(p) && !hist_matches(p) && !g_live.lattice.ok) cap = "Schrodinger";
    if (page == PAGE_PATHS || page == PAGE_DIAG) cap = "long Metropolis chain";
    if (page == PAGE_CLASS) cap = "outside the turning points";
    if (page == PAGE_POTS) cap = "ground energy";
    if (page == PAGE_HBAR || page == PAGE_CONV) cap = "analytic ground state";
    text(g_ui, cap, r.x + 14, r.y + 72, 14, COL_MUTED);
    if (r.height < 140) return;
    char line[96];
    const Spectrum& se_show = (live_num && g_live.se.ok) ? g_live.se : p.schrodinger;
    snprintf(line, sizeof line, "Schrodinger   %.5f", se_show.ok ? interp(se_show.x, se_show.psi2, 1.0) : 0.0);
    text(g_ui, line, r.x + 14, r.y + 100, 15, COL_INK);
    double analytic_show = (live_num && g_live.pot == 0) ? ho_ground_psi2(1.0, g_live.mass, g_live.hbar, 1.0) : p.analytic1;
    if (analytic_show > 0.0) snprintf(line, sizeof line, "analytic      %.5f", analytic_show);
    else snprintf(line, sizeof line, "analytic      n/a");
    text(g_ui, line, r.x + 14, r.y + 122, 15, COL_INK);
    const HoLattice& lat_show = live_num ? g_live.lattice : p.lattice;
    if (lat_show.ok) snprintf(line, sizeof line, "lattice       %.5f", gaussian_psi2(1.0, lat_show.alpha));
    else snprintf(line, sizeof line, "lattice       n/a");
    text(g_ui, line, r.x + 14, r.y + 144, 15, COL_INK);
    snprintf(line, sizeof line, "Metropolis    %.4f +/- %.4f", p.metro_small.psi2_at_1, p.metro_small.err_at_1);
    text(g_ui, line, r.x + 14, r.y + 166, 15, COL_INK);
}

int main(int argc, char** argv) {
    bool shots = false;
    bool histshot = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--shots") == 0) shots = true;
        if (strcmp(argv[i], "--hist") == 0) histshot = true;
    }

    FILE* log = shots ? fopen("build/runlog.txt", "w") : nullptr;
    if (log) {
        fprintf(log, "start shots=%d\n", shots ? 1 : 0);
        fflush(log);
    }
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT | FLAG_VSYNC_HINT | FLAG_WINDOW_HIGHDPI);
    InitWindow(1480, 880, "Feynman path integral");
    if (log) {
        fprintf(log, "window ready=%d size=%d x %d close=%d\n", IsWindowReady() ? 1 : 0, GetScreenWidth(), GetScreenHeight(), WindowShouldClose() ? 1 : 0);
        fflush(log);
    }
    SetWindowMinSize(1100, 700);
    SetWindowPosition(60, 30);
    SetTargetFPS(60);

    g_ui = LoadFontEx("C:/Windows/Fonts/segoeui.ttf", 48, nullptr, 0);
    g_code = LoadFontEx("C:/Windows/Fonts/consola.ttf", 40, nullptr, 0);
    if (g_ui.texture.id) SetTextureFilter(g_ui.texture, TEXTURE_FILTER_BILINEAR);
    if (g_code.texture.id) SetTextureFilter(g_code.texture, TEXTURE_FILTER_BILINEAR);

    refresh_exact();
    std::thread worker(worker_main);
    request_job(1, 7, 4.0, 1.0, 1.0, 0);

    int page = PAGE_MILE;
    int shot_i = -1;
    while (!WindowShouldClose()) {
        std::shared_ptr<Pack> pack;
        {
            std::lock_guard<std::mutex> lock(g_mu);
            pack = g_pack;
        }
        static int shot_wait = 0;
        if (shots && pack) {
            shot_wait++;
            if (GetScreenWidth() < 1200 && shot_wait < 45) {
                SetWindowSize(1480, 880);
                SetWindowPosition(40, 24);
            }
            if (shot_i < 0 && (GetScreenWidth() >= 1200 || shot_wait > 60)) shot_i = 0;
        }
        if (shots && shot_i >= 0 && shot_i < PAGE_COUNT) page = shot_i;

        int key = GetCharPressed();
        if (key >= '1' && key < '1' + PAGE_COUNT && !shots) {
            page = key - '1';
            g_scroll = 0;
        }
        if (!shots && (IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_KP_ENTER))) launch_run();
        if (g_play && !shots && pack && !g_busy.load()) {
            if (g_set_N >= 48) g_play = false;
            else if (sample_matches(*pack)) {
                g_set_N += 1;
                refresh_exact();
                kick_live();
            }
        }

        BeginDrawing();
        ClearBackground(COL_BG);
        int W = GetScreenWidth();
        int H = GetScreenHeight();
        DrawRectangle(0, 0, W, 4, COL_TEAL);
        text(g_ui, "FEYNMAN", 22, 16, 26, COL_INK);
        text(g_ui, "path integral engine", 168, 24, 16, COL_MUTED);
        if (g_busy) text(g_ui, "updating", (float)W - 110, 22, 16, COL_GOLD);

        if (!pack) {
            DrawRectangleRounded(Rectangle{W * 0.5f - 280, H * 0.5f - 70, 560, 140}, 0.08f, 8, COL_PANEL);
            text(g_ui, stage_name(g_stage), W * 0.5f - 250, H * 0.5f - 20, 22, COL_INK);
            float t = (float)(GetTime() * 0.8);
            DrawRectangle(W * 0.5f - 200, H * 0.5f + 28, 400, 4, COL_LINE);
            DrawRectangle(W * 0.5f - 200, H * 0.5f + 28, 80 + 40 * (1 + sinf(t * 3)), 4, COL_TEAL);
        } else {
        float top = 56;
        float bot = 30;
        float gap = 14;
        float sideW = 168.0f;
        float rightW = 276.0f;
        float barH = 112.0f;
        float bodyH = (float)H - top - bot;
        float centerW = (float)W - 32.0f - sideW - rightW - 2.0f * gap;
        if (centerW < 520.0f) {
            rightW = 0;
            centerW = (float)W - 32.0f - sideW - gap;
        }
        if (bodyH < barH + 280.0f) barH = std::max(96.0f, bodyH - 280.0f);
        float plotH = bodyH - barH - gap;
        Rectangle side{16, top, sideW, bodyH};
        card(side);
        float step = std::min(40.0f, (side.height - 48.0f) / (float)PAGE_COUNT);
        for (int i = 0; i < PAGE_COUNT; ++i) {
            Rectangle b{side.x + 8, side.y + 10 + i * step, sideW - 16, step - 6};
            if (button(b, kPages[i], page == i) && !shots) {
                page = i;
                g_scroll = 0;
            }
        }
        Rectangle plot{16 + sideW + gap, top, centerW, plotH};
        Rectangle bar{plot.x, plot.y + plotH + gap, centerW, barH};
        if (rightW > 0) {
            Rectangle right{plot.x + centerW + gap, top, rightW, bodyH};
            draw_numbers(Rectangle{right.x, right.y, right.width, 188}, page, *pack);
            Rectangle script{right.x, right.y + 202, right.width, right.height - 202};
            draw_script(script, script_for(page));
        }
        draw_controls(bar);
        BeginScissorMode((int)plot.x, (int)plot.y, (int)plot.width, (int)plot.height);
        if (page == PAGE_MILE) draw_milestone(plot, *pack);
        else if (page == PAGE_CONV) draw_convergence(plot, *pack);
        else if (page == PAGE_PATHS) draw_paths(plot, *pack);
        else if (page == PAGE_DIAG) draw_diagnostics(plot, *pack);
        else if (page == PAGE_GAP) draw_gap(plot, *pack);
        else if (page == PAGE_CLASS) draw_classical(plot, *pack);
        else if (page == PAGE_HBAR) draw_hbar(plot, *pack);
        else if (page == PAGE_POTS) draw_pots(plot, *pack);
        else if (page == PAGE_RESEARCH) draw_research(plot, *pack);
        else if (page == PAGE_LIBRARY) draw_library(plot, *pack);
        else draw_design(plot, *pack);
        EndScissorMode();
    }
    DrawRectangle(0, H - 26, W, 26, Color{12, 14, 20, 255});
    text(g_ui, "N goes up to a million. To 48 the curve follows. Above that, Run writes data/histogram.csv. Keys 1-9 change page.", 16, (float)H - 20, 14, COL_MUTED);

        EndDrawing();
        static int shot_phase = 0;
        if (shots && pack && shot_i >= 0 && shot_i < PAGE_COUNT) {
            if (shot_phase == 1) {
                char name[64];
                snprintf(name, sizeof name, "shot_%d.png", shot_i);
                TakeScreenshot(name);
                shot_i++;
                shot_phase = 0;
            } else {
                shot_phase = 1;
            }
        }
        if (shots && shot_i >= PAGE_COUNT) break;
        if (histshot && pack) {
            static int hist_phase = 0;
            if (GetScreenWidth() < 1200) {
                SetWindowSize(1480, 880);
                SetWindowPosition(40, 24);
            }
            if (hist_phase == 0 && !g_busy.load()) {
                g_set_N = 1000000;
                g_set_T = 4;
                g_set_hbar = 1;
                g_set_mass = 1;
                g_set_pot = 0;
                refresh_exact();
                request_job(4, 1000000, 4.0, 1.0, 1.0, 0);
                hist_phase = 1;
            } else if (hist_phase == 1 && !g_busy.load() && pack->hist.ok && pack->hist.N == 1000000) {
                page = PAGE_MILE;
                hist_phase = 2;
            } else if (hist_phase == 2) {
                hist_phase = 3;
            } else if (hist_phase == 3) {
                TakeScreenshot("shot_hist.png");
                break;
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(g_mu);
        g_quit = true;
    }
    g_cv.notify_one();
    worker.join();
    if (g_world.id) UnloadTexture(g_world);
    if (g_ui.texture.id) UnloadFont(g_ui);
    if (g_code.texture.id) UnloadFont(g_code);
    CloseWindow();
    if (log) {
        fprintf(log, "exit\n");
        fclose(log);
    }
    return 0;
}
