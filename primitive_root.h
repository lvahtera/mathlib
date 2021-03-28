#pragma once
#include "primes_t.h"
#include "mod_a_t.h"

template <typename T>
T primitive_root(T p, bool prime = true);

template <typename T>
inline bool pr_subf(T k, T p, const std::vector<T>& factors);

/* Primitive root modulo n
   returns 0 if it doesn't exist
   */
template <typename T>
T primitive_root(T n, bool prime) {
    T phi = prime ? n - 1 : eulers_totient(n);
    std::vector<T> factors = factorise(phi);
    for (T k = 2; k <= n; ++k) {
        if (pr_subf(k, n, factors)) return k;
    }
    return -1;
}

template <typename T>
inline bool pr_subf(T k, T n, const std::vector<T>& factors) {
    T prev = 0;
    for (auto f: factors) {
        if (f == prev) continue;
        prev = f;
        if (mod_exp(k, (n - 1) / f, n) == 1) return false;
    }
    return true;
}
