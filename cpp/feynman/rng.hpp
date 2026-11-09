#pragma once
#include <cstdint>
#include <cmath>

namespace feynman {

inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

struct Rng {
    uint64_t s0, s1, s2, s3;
    double spare;
    bool has_spare;

    explicit Rng(uint64_t seed) : spare(0.0), has_spare(false) {
        s0 = splitmix(seed);
        s1 = splitmix(seed);
        s2 = splitmix(seed);
        s3 = splitmix(seed);
    }

    static uint64_t splitmix(uint64_t& x) {
        x += 0x9e3779b97f4a7c15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    uint64_t next_u64() {
        uint64_t result = rotl64(s1 * 5, 7) * 9;
        uint64_t t = s1 << 17;
        s2 ^= s0;
        s3 ^= s1;
        s1 ^= s2;
        s0 ^= s3;
        s2 ^= t;
        s3 = rotl64(s3, 45);
        return result;
    }

    double uniform() {
        return (next_u64() >> 11) * (1.0 / 9007199254740992.0);
    }

    double normal() {
        if (has_spare) {
            has_spare = false;
            return spare;
        }
        double u1 = uniform();
        double u2 = uniform();
        if (u1 < 1e-300) u1 = 1e-300;
        double r = std::sqrt(-2.0 * std::log(u1));
        double theta = 6.28318530717958647692 * u2;
        spare = r * std::sin(theta);
        has_spare = true;
        return r * std::cos(theta);
    }
};

}
