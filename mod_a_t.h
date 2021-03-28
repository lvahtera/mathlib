#pragma once
#include <cstdint>
#include <algorithm>
#include <limits>
#include <vector>

using u32 = std::uint32_t;

template <typename T>
T mod_mult(T a, T b, T mod);

template <typename T1, typename T2>
T1 mod_exp(T1 a, T2 e, T1 mod);

/* Modular multiplication
   0 <= a, b < 2^64, 0 < mod < 2^63 (for T_t)
   */
template <typename T>
T mod_mult(T a, T b, const T mod) {
    if (a >= mod || a < 0) a = (a % mod + mod) % mod;
    if (b >= mod || b < 0) b = (b % mod + mod) % mod;
    if ((a | b) < std::numeric_limits<u32>::max()) {
        return (a * b) % mod;
    }

    T r = 0;         // remainder
    T m = mod >> 1;  // for n < m, if 2 * n > mod -> 2n % mod

    if (b > a) std::swap(a, b);

    while (b) {
        // multiplies "a" by "b" bit by bit, cf long multiplication
        if (b & 1) {
            r += a;
            if (r >= mod) r -= mod;
        }
        b >>= 1;
        a = a > m ? (a << 1) - mod : a << 1;
    }
    return r;
}

/* Modular exponentiation
   Exponentiation by squaring
   0 <= a, e < 2^64, 0 < mod < 2^63 (for uint64_t)
   */
template <typename T1, typename T2>
T1 mod_exp(T1 a, T2 e, const T1 mod) {
    T1 r = 1;
    if (a >= mod || a < 0) a = (a % mod + mod) % mod;

    if (mod < std::numeric_limits<u32>::max()) {
        while (e) {
            if (e & 1) r = r * a % mod;
            e >>= 1;
            a = a * a % mod;
        }
    }
    else {
        while (e) {
            if (e & 1) r = mod_mult(a, r, mod);
            e >>= 1;
            a = mod_mult(a, a, mod);
        }
    }
    return r;
}
