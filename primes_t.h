#pragma once
#include "mod_a_t.h"
#include "misc_al_t.h"
#include <iterator>
#include <algorithm>
#include <vector>

using u64 = std::uint64_t;
using u8 = std::uint8_t;

template <typename T>
std::vector<T> gen_primes(T lim);

template <typename T>
std::vector<T> segmented_sieve(T lim);

template <typename T>
std::vector<T> wheel_factor(const std::vector<T>& primes, T limit, const std::vector<T>&& cur_wheel = std::vector<T>());

template <typename T>
std::vector<T> wheel_rotation(int rotations, T wheel_d, int prime_c, const std::vector<T>&& cur_wheel);

template <typename T>
std::vector<T> wheel_incr(const std::vector<T>& primes);

template <typename T>
std::vector<T> sieve_of_eratos(const std::vector<T>& input, int basis_size = 0, int lim = 0);

template <typename T>
std::vector<T> factorise(T n);

template <typename T>
bool is_prime(T n);

template <typename T>
bool trial_div(T n);

template <typename T>
bool miller_rabin(T n);

template <typename T>
bool mr_is_composite(T n, T d, int s, T a);

template <typename T>
T eulers_totient(T n);

/* Prime generator
   Uses the functions below to generate all primes below T lim
   */
template <typename T>
std::vector<T> gen_primes(T lim) {
    return segmented_sieve(lim);
}

/* Segmented sieve of Eratosthenes with a small wheel.
 * Based on Kim Walisch's code at
 * https://github.com/kimwalisch/primesieve/wiki/Segmented-sieve-of-Eratosthenes
 */

template <typename T>
std::vector<T> segmented_sieve(T lim) {
    if (lim < 2) return std::vector<T>{};
    if (lim < 3) return std::vector<T>{2};
    if (lim < 5) return std::vector<T>{2, 3};

    constexpr int segment_size = 524288; // L2 cache size
    T sqrt_lim = isqrt(lim);

    std::vector<u8> is_prime(sqrt_lim + 1, true);
    std::vector<u8> sieve(segment_size);
    std::vector<T> sieve_primes = {};
    std::vector<T> multiples = {};

    std::vector<T> primes(T(1.26 * (lim / log(lim))) + 1);
    primes[0] = 2; primes[1] = 3; primes[2] = 5;

    T a = 3;
    T n = 7;
    T pc = 3; // # of primes

    // 2*3*5 sieving wheel
    u8 o = 0;
    constexpr u8 offset[] = {4, 2, 4, 2, 4, 6, 2, 6};

    for (T low = 0; low <= lim; low += segment_size) {
        std::fill(sieve.begin(), sieve.end(), true);
        T high = low + segment_size - 1;
        high = std::min(high, lim);

        while (a * a <= high) {
            if (is_prime[a]) {
                sieve_primes.emplace_back(a);
                multiples.emplace_back(a * a - low);
                for (T b = a * a; b <= sqrt_lim; b += a) {
                    is_prime[b] = false;
                }
            }
            a += 2;
        }

        for (int i = 0; i < sieve_primes.size(); ++i) {
            T j = multiples[i];
            for (T k = 2 * sieve_primes[i]; j < segment_size; j += k) sieve[j] = false;
            multiples[i] = j - segment_size;
        }

        while (n <= high) {
            if (sieve[n - low]) primes[pc++] = n;
            n += offset[o++];
            o &= 0x7;
        }
    }
    if (primes.size() > pc) primes.resize(pc);
    return primes;
}

/* Wheel Factorization.
   Given a basis (a few prime numbers) generates a list of integers that are
   coprime (and mostly prime) with all the numbers of the basis.
*/
template <typename T>
std::vector<T> wheel_factor(const std::vector<T>& primes, T limit,
                            const std::vector<T>&& cur_wheel) {
    if (!limit) {
        return primes;
    }
    T wheel_d{1};
    for (auto p: primes) wheel_d *= p;
    T rotations = limit / wheel_d + 1;

    if (!cur_wheel.empty()) {
        return wheel_rotation(rotations, wheel_d, primes.size(), std::move(cur_wheel));
    }
    std::vector<T> wheel(primes);
    wheel.reserve(wheel_d);
    // Sieve out the composites, only guaranteed to be coprime to the basis
    for (T i = primes.back() + 1; i < wheel_d; ++i) {
        bool is_composite = false;
        for (auto prime: primes) {
            if (!(i % prime)) {
                is_composite = true;
                break;
            }
        }
        if (!is_composite) {
            wheel.push_back(i);
        }
    }
    rotations--;
    if (!rotations) {
        return wheel;
    }
    else return wheel_rotation(rotations, wheel_d, primes.size(), std::move(wheel));
}

// additional layers
template <typename T>
std::vector<T> wheel_rotation(int rotations, T wheel_d, int prime_c,
                                const std::vector<T>&& cur_wheel) {
    std::vector<T> wheel(std::move(cur_wheel));
    int spoke_count = std::distance(wheel.begin(), std::upper_bound(wheel.begin(), wheel.end(), wheel_d));
    spoke_count -= prime_c - 1;

    /* Numbers along the spokes are coprime to the basis, i.e.
    wheel_d + wheel_spokes[0] is relatively prime to the basis etc
    */
    std::vector<T> wheel_spokes(spoke_count);
    wheel_spokes[0] = 1;
    for (int i = 1; i < spoke_count; ++i) {
        wheel_spokes[i] = wheel[prime_c + i - 1];
    }
    T n = (wheel.back() / wheel_d) + 1;
    wheel.reserve(wheel.size() + n * wheel_spokes.size());

    while (rotations) {
        for (auto relative_prime: wheel_spokes) {
            wheel.emplace_back(n * wheel_d + relative_prime);
        }
        ++n;
        rotations--;
    }
    return wheel;
}

/* Use this to generate potential primes if you don't want to create a
   large array. See "spokes" above.
   */
template <typename T>
std::vector<T> wheel_incr(const std::vector<T>& primes) {
    T wheel_d{1};
    for (auto p: primes) wheel_d *= p;
    std::vector<T> increments{};
    std::vector<T> wheel = wheel_factor(primes, 2 * wheel_d);
    for (int i = primes.size(); wheel[i] <= wheel_d + 1; ++i) {
        increments.emplace_back(wheel[i+1] - wheel[i]);
    }
    return increments;
}

/* Sieve of Eratosthenes
   No bells and whistles
   basis_size = # of numbers that are 100% prime
 */
template <typename T>
std::vector<T> sieve_of_eratos(const std::vector<T>& nr_list, int basis_size, int lim) {
    if (lim == 0) lim = nr_list.back();
    std::vector<bool> candidates(lim + 1, false);
    for (auto p: nr_list) if (p <= lim) candidates[p] = true;
    T multiples = 0; // for prime p: p^2, p^2 + p, p^2 + 2p...
    for (T i = nr_list[basis_size]; (multiples = i * i) < lim; i = nr_list[++basis_size]) {
        if (candidates[i]) {
            while (multiples <= lim) {
                candidates[multiples] = false;
                multiples += i;
            }
        }
    }
    std::vector<T> primes{};
    for (T i = 0; i < candidates.size(); ++i) {
        if (candidates[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

/* Factorises given number.
   12 -> {2, 2, 3}
   */
template <typename T>
std::vector<T> factorise(T n) {
    if (n == 1) return {};
    std::vector<T> factors{};
    static const std::vector<T> primes = {2, 3, 5, 7};
    static const std::vector<T> increments = wheel_incr(primes);
    for (auto p: primes) {
        while (n % p == 0) {
            factors.push_back(p);
            n /= p;
        }
    }

    int i{0};
    for (T div = 11; div * div <= n; div += increments[i++]) {
        while (n % div == 0) {
            factors.push_back(div);
            n /= div;
        }
        if (i == increments.size()) {
            i = 0;
        }
    }
    // only true if n is prime
    if (n > 1) {
        factors.push_back(n);
    }

    return factors;
}

template <typename T>
bool is_prime(T n) {
    if (n < 2) {
        return false;
    }
    if (n < 321503171) {
        return trial_div(n);
    }
    if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0 ||
        n % 11 == 0 || n % 13 == 0 || n % 17 == 0 || n % 19 == 0) {
        return false;
    }
    return miller_rabin(n);
}

/* Simple primality test
   */
template <typename T>
bool trial_div(T n) {
    if (n < 2) {
        return false;
    }
    static const std::vector<T> primes = {2, 3, 5, 7, 11};
    static const std::vector<T> increments = wheel_incr(primes);
    for (auto p: primes) {
        if (n == p) {
            return true;
        }
        else if (n % p == 0) {
            return false;
        }
    }

    int i{0};
    for (T div = 13; div * div <= n; div += increments[i++]) {
        if (n % div == 0) {
            return false;
        }
        if (i == increments.size()) {
            i = 0;
        }
    }
    return true;
}

/* Miller-Rabin
   Works if n < 2^64
   */
template <typename T>
bool miller_rabin(T n) {
    if (!(n & 1) || n < 2) {
        return false;
    }
    T d = n - 1;
    int r = 0;
    while (!(d & 1)) {
        ++r;
        d >>= 1;
    }

    std::vector<T> arr = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

    for (auto a: arr) {
        if (n == a) {
            return true;
        }
        else if (mr_is_composite(n, d, r, a)) {
            return false;
        }
    }

    return true;
}

template <typename T>
bool mr_is_composite(T n, T d, int r, T a) {
    T x = mod_exp(a, d, n);
    if (x == 1 || x == n - 1) {
        return false;
    }
    for (int i = 1; i < r; ++i) {
        x = mod_exp<T>(x, 2, n);
        if (x == n - 1) {
            return false;
        }
    }

    return true;
}

template <typename T>
T eulers_totient(T n) {
    std::vector<T> fact = factorise(n);
    T num = 1;
    T den = 1;
    T prev = 0;

    for (auto& p: fact) {
        if (p = prev) continue;
        num *= p - 1;
        den *= p;
        prev = p;
    }
    return n * num / den;
}
