#pragma once
#include <string>
#include <cmath>

namespace feynman {

enum PotKind { HO = 0, ANHARM = 1, DOUBLEWELL = 2, MORSE = 3, QUARTIC = 4 };

struct Potential {
    int kind;
    double a;
    double b;
    double c;

    double operator()(double x) const {
        switch (kind) {
        case HO:
            return 0.5 * a * x * x;
        case ANHARM:
            return 0.5 * a * x * x + b * x * x * x * x;
        case DOUBLEWELL: {
            double d = x * x - b * b;
            return a * d * d;
        }
        case MORSE: {
            double e = 1.0 - std::exp(-b * (x - c));
            return a * e * e;
        }
        case QUARTIC:
            return a * x * x * x * x;
        default:
            return 0.5 * x * x;
        }
    }

    const char* name() const {
        switch (kind) {
        case HO: return "ho";
        case ANHARM: return "anharm";
        case DOUBLEWELL: return "doublewell";
        case MORSE: return "morse";
        case QUARTIC: return "quartic";
        default: return "ho";
        }
    }

    const char* label() const {
        switch (kind) {
        case HO: return "harmonic oscillator";
        case ANHARM: return "anharmonic oscillator";
        case DOUBLEWELL: return "double well";
        case MORSE: return "Morse";
        case QUARTIC: return "quartic";
        default: return "harmonic oscillator";
        }
    }
};

inline Potential parse_potential(const std::string& s) {
    if (s == "anharm" || s == "anharmonic") return Potential{ANHARM, 1.0, 0.25, 0.0};
    if (s == "doublewell" || s == "dw") return Potential{DOUBLEWELL, 0.8, 1.5, 0.0};
    if (s == "morse") return Potential{MORSE, 3.0, 0.65, 0.0};
    if (s == "quartic") return Potential{QUARTIC, 0.25, 0.0, 0.0};
    return Potential{HO, 1.0, 0.0, 0.0};
}

inline double ho_omega(const Potential& V, double mass) {
    if (V.kind != HO || mass <= 0.0) return 1.0;
    return std::sqrt(V.a / mass);
}

inline double ho_ground_psi2(double x, double mass, double hbar, double k) {
    double omega = std::sqrt(k / mass);
    double alpha = mass * omega / hbar;
    return std::sqrt(alpha / 3.14159265358979323846) * std::exp(-alpha * x * x);
}

inline double ho_finite_T_psi2(double x, double mass, double hbar, double k, double T) {
    double omega = std::sqrt(k / mass);
    double th = std::tanh(0.5 * omega * T);
    double alpha = (mass * omega / hbar) * th;
    return std::sqrt(alpha / 3.14159265358979323846) * std::exp(-alpha * x * x);
}

}
