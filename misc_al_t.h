#pragma once
#include <cstdint>
#include <cmath>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;

template <typename T>
T binpow_fibonacci(int n);

template <typename T>
T binpow_narayana(int n);

template <typename T>
T fast_fibonacci(u64 n);

template <typename T1, typename T2>
T1 bpow(T1 a, T2 exp);

template <typename T>
T xgcd(T a, T b, i64& x, i64& y);

template <typename T>
T crt(T a, T m, T b, T n);

template <typename T1, typename T2>
T1 nck(T1 n, T2 k);

template <typename T1, typename T2>
T1 iroot(T1 x, T2 n) {
    T1 res = T1(pow(x, 1.0 / n));
    if (bpow(res, n) > x) --res;
    else if (bpow(res + 1, n) <= x) ++res;
    return res;
}

template <typename T>
T isqrt(T n) { return iroot(n, 2); }

template <typename T>
T icbrt(T n) { return iroot(n, 3); }

// index of the MSB
inline unsigned int msb(u64 n) {
    return n ? 63 - __builtin_clzll(n) : 0;
}

/* Utilises matrices to rapidly calculate F_n. Starts from the bit after MSB.
   */
template <typename T>
T binpow_fibonacci(int n) {
    if (n < 2) {
        return n == 1 ? 1 : 0;
    }
    T a = 1, b = 1, c = 0; // matrix {{1, 1},{1,0}} = {{a, b}, {b, c}}
    int i = 1 << msb(n);
    while (i >>= 1) {
        T t = b * b;
        b = b * (a + c);
        a = a * a + t;
        c = t + c * c;
        if (n & i) {
            c = b;
            b = a;
            a = a + c; // actually a + b
        }
    }
    return b;
}

/* Rapidly counts Narayana's cows = N_n = N_n-1 + N_n-3
   */
template <typename T>
T binpow_narayana(int n) {
    // {{1,0,1},{1,0,0},{0,1,0}} = {{a,c,b},{b,d,c},{c,e,d}}
    T a = 1, b = 1, c = 0;
    int i = 1 << msb(n); // the MSB, safe because n > 1
    while (i >>= 1) {
        T t1 = 2 * b * c;
        T t2 = b * b;
        b = c * c + 2 * a * b - t2;
        c = t2 + 2 * a * c - t1;
        a = a * a + t1;
        if (n & i) {
            T t_c = c;
            c = b;
            b = a;
            a = a + t_c;
        }
    }
    return a;
}

/* See: https://en.wikipedia.org/wiki/Fibonacci_number#Matrix_form
   */
template <typename T>
T fast_fibonacci(u64 n) {
    if (n < 2) {
        return T(n);
    }
    T fn = 1;
    T fnm = 0;
    int i = 1 << msb(n);
    while (i >>= 1) {
        T f2nm = fn * fn + fnm * fnm;
        T f2n = fn * (2 * fnm + fn);
        fn = n & i ? f2n + f2nm : f2n;
        fnm = n & i ? f2n : f2nm;
    }
    return fn;
}

template <typename T1, typename T2>
T1 bpow(T1 a, T2 exp) {
    T1 b = 1;
    while (exp > 1) {
        if (exp & 1) {
            b *= a;
        }
        a *= a;
        exp >>= 1;
    }
    return a * b;
}

template <typename T>
T xgcd(T a, T b, i64& x, i64& y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }

    i64 x1, y1;
    T d = xgcd(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

template <typename T>
T crt(T a, T m, T b, T n) {
    i64 x, y;
    xgcd(m, n, x, y);
    T mod = m * n;
    return (mod + (a + (b - a) * x * m) % mod) % mod;
}

template <typename T1, typename T2>
T1 nck(T1 n, T2 k) {
    if (n < k) return 0;
    if (n == k || k == 0) return 1;
    u64 res = k >= n - k ? k + 1 : n - k + 1;
    for (u64 i = res + 1, j = 2; i <= n; ++i, ++j) {
        res *= i;
        res /= j;
    }
    return res;
}
