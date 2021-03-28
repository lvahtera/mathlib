#pragma once
#include "mod_a_t.h"
#include <cstdint>
#include <iostream>

/* Galois Field of order Mod, which (currently) has to be a prime.
 * If Mod is not a prime number, modular inverse might not be defined
 * and division might (will) give inaccurate results.
 * Other operations still guaranteed to work in Z/nZ.
 * */

using u64 = std::uint64_t;
using i64 = std::int64_t;

template <i64 Mod> class Field {
// TODO rhs % mod == 0 -> no modinverse
private:
    i64 a;

public:
    constexpr Field(const i64 n = 0) noexcept : a(n % Mod) {}
    constexpr Field(const Field& n) noexcept : a(n.a) {}

    constexpr i64 value() const noexcept { return (a + Mod) % Mod; }
    explicit constexpr operator i64() const noexcept { return a; }
    explicit constexpr operator int() const noexcept { return a; }
    explicit constexpr operator bool() const noexcept { return a; }

    constexpr Field pow(i64 exp) const noexcept { return Field(mod_exp(a, exp, Mod)); }

    constexpr Field& operator+=(const Field& rhs) noexcept {
        a += rhs.a;
        if (a >= Mod) a -= Mod;
        return *this;
    }

    constexpr Field& operator-=(const Field& rhs) noexcept {
        if (a < rhs.a) a += Mod;
        a -= rhs.a;
        return *this;
    }

    constexpr Field& operator*=(const Field& rhs) noexcept {
        a = mod_mult<i64>(a, rhs.a, Mod);
        return *this;
    }

    constexpr Field& operator/=(const Field& rhs) noexcept {
        *this *= mod_exp(rhs.a, Mod - 2, Mod);
        return *this;
    }

    constexpr Field& operator+=(i64 rhs) noexcept {
        if (Mod <= rhs || rhs < 0) rhs = ((rhs % Mod) + Mod) % Mod;
        a += rhs;
        if (a >= Mod) a -= Mod;
        return *this;
    }

    constexpr Field& operator-=(i64 rhs) noexcept {
        if (Mod <= rhs || rhs < 0) rhs = ((rhs % Mod) + Mod) % Mod;
        if (a < rhs) a += Mod;
        a -= rhs;
        return *this;
    }

    constexpr Field& operator*=(i64 rhs) noexcept {
        if (Mod <= rhs || rhs < 0) rhs = ((rhs % Mod) + Mod) % Mod;
        a = mod_mult<i64>(a, rhs, Mod);
        return *this;
    }

    constexpr Field& operator/=(i64 rhs) noexcept {
        if (Mod <= rhs || rhs < 0) rhs = ((rhs % Mod) + Mod) % Mod;
        *this *= mod_exp<i64>(rhs, Mod - 2, Mod);
        return *this;
    }

    constexpr Field operator+(const Field& rhs) const noexcept { return Field(*this) += rhs; }
    constexpr Field operator-(const Field& rhs) const noexcept { return Field(*this) -= rhs; }
    constexpr Field operator*(const Field& rhs) const noexcept { return Field(*this) *= rhs; }
    constexpr Field operator/(const Field& rhs) const noexcept { return Field(*this) /= rhs; }

    constexpr Field operator-() const noexcept { return Field(Mod - a); }
    constexpr bool operator==(const Field& rhs) const noexcept { return a == rhs.a; }
    constexpr bool operator!=(const Field& rhs) const noexcept { return a != rhs.a; }

    constexpr Field operator+(const i64 rhs) const noexcept { return Field(*this) += rhs; }
    constexpr Field operator-(const i64 rhs) const noexcept { return Field(*this) -= rhs; }
    constexpr Field operator*(const i64 rhs) const noexcept { return Field(*this) *= rhs; }
    constexpr Field operator/(const i64 rhs) const noexcept { return Field(*this) /= rhs; }
};

template <i64 Mod>
std::istream& operator>>(std::istream& in, Field<Mod>& f) {
    i64 x;
    in >> x;
    f = Field<Mod>(x);
    return in;
}

template <i64 Mod>
std::ostream& operator<<(std::ostream& out, const Field<Mod>& f) {
    out << f.value();
    return out;
}

template <i64 Mod>
Field<Mod> pow(const Field<Mod>& base, i64 exp) {
    return Field<Mod>(mod_exp(base, exp, Mod));
}
