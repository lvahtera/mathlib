#pragma once
#include "misc_al_t.h"
#include "primes_t.h"
#include <algorithm>
#include <cstdint>
#include <vector>
#include <cmath>
#include <array>
#include <unordered_map>

using u64 = std::uint64_t;

template <typename T1, typename T2>
u64 pcount_phi(const std::vector<T1>& primes, T2 n, u64 a);

template <typename T1, typename T2>
u64 lehmer_pi(const std::vector<T1>& primes, T2 n);

template <typename T1, typename T2>
u64 PX(const std::vector<T1>& primes, T2 n, u64 a = 0, u64 x = 2);

template <typename T>
u64 primecount(T n);

template <typename T>
u64 primecount(T n) {
    if (n < 2) return 0;
    constexpr u64 pcount = 1000*1000*100;
    static const std::vector<T> primes = gen_primes<u64>(pcount);
    return lehmer_pi(primes, n);
}

/* phi(n, a) returns the number of primes below n minus the first a primes.
 * */
template <typename T1, typename T2>
u64 pcount_phi(const std::vector<T1>& primes, T2 n, u64 a) {
    // cache for certain values of a and x that are commonly visited
    constexpr size_t max_a = 0xff;
    constexpr u64 max_n = 0xffff;
    static std::array<std::array<u64, max_n>, max_a> phi_cache = {};

    if (n < 1) return 0;
    else if (n < a) return 1;
    else if (a < 1) return n;
    else if (a < max_a && n < max_n && phi_cache[a][n]) return phi_cache[a][n];

    u64 phi_sum = pcount_phi(primes, n, a - 1) - pcount_phi(primes, n / primes[a-1], a - 1);
    if (n < max_n && a < max_a) phi_cache[a][n] = phi_sum;
    return phi_sum;
}

/*  Meissel-Lehner prime counting function
 *  can be easily transformed into Meissel's or Legendre's original function(s)
 *  as they're all based on each other.
 *  Utilises caches to store commonly visited values, one for small pi(n) to
 *  avoid traversing the vector every time, and the other to avoid calculating
 *  pi(n) multiple times for large n. lpc is only useful if this function is
 *  called dozens or hundreds of times with large n. Enable/disable as needed.
 *  */
template <typename T1, typename T2>
u64 lehmer_pi(const std::vector<T1>& primes, T2 n) {
    constexpr u64 max_n = 0xffff;
    const u64 large_pi = primes.back();
    static std::unordered_map<u64, u64> lpc = {};   // cache for large n
    static std::array<u64, max_n> pi_cache = {};    // cache for n < 0xffff

    if (n > large_pi && lpc.count(n)) return lpc[n];
    if (n < max_n && pi_cache[n]) return pi_cache[n];
    if (n <= primes.back()) {
        u64 prime_c = std::upper_bound(primes.begin(), primes.end(), n) - primes.begin();
        if (n < max_n) pi_cache[n] = prime_c;
        return prime_c;
    }

    u64 root = iroot(n, 4);
    u64 a = lehmer_pi(primes, root);
    u64 pi = pcount_phi(primes, n, a) + a - 1 - PX(primes, n, a, 2) - PX(primes, n, a, 3);

    if (n > large_pi) lpc[n] = pi;

    return pi;
}

/* PX(n, a) calculates the number of k-almost primes between the p_ath prime
 * and n.
 * */
template <typename T1, typename T2>
u64 PX(const std::vector<T1>& primes, T2 n, u64 a, u64 x) {
    if (x == 0) return 1;
    if (x == 1) return lehmer_pi(primes, n);
    u64 sum = 0;
    u64 b = lehmer_pi(primes, iroot(n, x));
    if (x == 2) {
        for (u64 i = a; i < b; ++i) {
            sum += lehmer_pi(primes, n / primes[i]) - i;
        }
    }
    else {
        for (u64 i = a; i < b; ++i) {
            sum += PX(primes, n / primes[i], i, x - 1);
        }
    }

    return sum;
}
