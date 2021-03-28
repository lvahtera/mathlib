#pragma once
#include "mod_a_t.h"
#include <cstdint>

using uint64 = std::uint64_t;
using int64 = std::int64_t;

template <typename T>
inline uint64 legendre(T a, T p);

template <typename T>
int64 tonelli(T n, T p);

template <typename T>
inline uint64 legendre(T a, T p) {
    return mod_exp(a, (p - 1) / 2, p);
}

/* Tonelli-Shanks algorithm 
   See wiki
   */
template <typename T>
int64 tonelli(T n, T p) {
    if (legendre(n, p) != 1) {
        return -1;
    }
    T Q = p - 1;
    T S = 0;
    while (!(Q & 1)) {
        Q /= 2;
        ++S;
    }

    T z = 2; // z is any quadratic nonresidue
    while (p - 1 != legendre(z, p) && z < p) {
        ++z;
    }

    T M = S;
    T c = mod_exp(z, Q, p);
    T t = mod_exp(n, Q, p);
    T R = mod_exp(n, (Q + 1) / 2, p);

    while (true) { 
        if (t == 0) return 0;
        else if (t == 1) return R;

        T t2 = (t * t) % p;
        int i = 1;
        while (t2 % p != 1 && i < M) {
            ++i;
            t2 = (t2 * t2) % p;
        }

        T b = mod_exp(c, 1ULL << (M - i - 1), p);
        R = (R * b) % p;
        c = (b * b) % p;
        t = (t * c) % p;
        M = i;
    }
    return 0;
}
