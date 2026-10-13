// Детерминированный ГСЧ: PCG32 на атом (потом переедет в compute-шейдер без изменений),
// splitmix64 для первичного засева.
#pragma once
#include <cmath>
#include <cstdint>

namespace decay::rng {

// PCG XSH-RR 64/32
struct Pcg32 {
    static constexpr uint64_t INC = 0xda3e39cb94b95bdbull; // нечётная константа потока

    static inline uint32_t next(uint64_t& s) {
        uint64_t old = s;
        s = old * 6364136223846793005ull + INC;
        uint32_t xs   = uint32_t(((old >> 18u) ^ old) >> 27u);
        uint32_t rot  = uint32_t(old >> 59u);
        return (xs >> rot) | (xs << ((0u - rot) & 31u)); // циклический сдвиг вправо на rot
    }

    static inline double u01(uint64_t& s) { // [0, 1), все 32 бита розыгрыша
        return double(next(s)) * (1.0 / 4294967296.0);
    }

    // время до распада ~ Exp(λ)
    static inline double expo(uint64_t& s, double lambda) {
        double u = u01(s);
        if (u <= 0.0) u = 1.0 / 9007199254740992.0; // защита от τ=0
        return -std::log(1.0 - u) / lambda;
    }
};

// первичный засев: каждое вызово даёт независимый 64-битный поток
inline uint64_t splitmix64(uint64_t& x) {
    x += 0x9e3779b97f4a7c15ull;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

} // namespace decay::rng