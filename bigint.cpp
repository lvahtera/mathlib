#include "bigint.h"
#include <iostream>
#include <vector>
#include <cassert>

/* constructors */
// default constructor
bigint::bigint(int n) {
    // -n when -min_value
    sign = n >= 0 ? POSITIVE : NEGATIVE;
    n = sign == POSITIVE ? n : -n;
    convert(n);
}

bigint::bigint(int64 n) {
    sign = n >= 0 ? POSITIVE : NEGATIVE;
    n = sign == POSITIVE ? n : -n;
    convert(n);
}

bigint::bigint(uint64 n) {
    sign = POSITIVE;
    convert(n);
}

// copy constructor
bigint::bigint(const bigint& n) {
    sign = n.sign;
    value = n.value;
}

/* assignment */
// copy assignment
bigint& bigint::operator=(const bigint& n) { // = default
    sign = n.sign;
    value = n.value;
    return *this;
}

// move assignment
bigint& bigint::operator=(bigint&& n) {
    sign = n.sign;
    value = std::move(n.value);
    return *this;
}

/* operators */
bigint& bigint::operator+=(const bigint& n) { // = default
    if (sign != n.sign) {
        return *this -= -n;
    
    }
    if (n.value.size() > value.size()) {
        value.resize(n.value.size(), 0);
    }

    uint64 carry = 0;
    for (int i = 0; i < value.size(); ++i) {
        uint64 sum = value[i] + carry + (n.value.size() > i ? n.value[i] : 0);
        carry = sum >> b_exp;   // faster than sum >= base;
        value[i] = sum & mask;  // faster than if (carry) sum -= base;
    
    }
    if (carry) {
        value.emplace_back(1);
    
    }
    return *this;
}

bigint& bigint::operator-=(const bigint& n) {
    if (sign != n.sign) {
        return *this += -n; // a, b >= 0; -a - b = -(a+b), a - (-b) = a+b
    }

    // a, b >= 0, a >= b; a-b
    int64 borrow = 0;
    if (gt_abs(*this, n)) {
        for (int i = 0; i < value.size(); ++i) {
            int64 diff = value[i] - (borrow + (n.value.size() > i ? n.value[i] : 0));
            borrow = diff < 0;
            value[i] = diff & mask;         // faster than if (borrow) diff += base;
        }

        // abs(*this) > abs(n) guarantees that value.empty() == false
        while (value.back() == 0) { 
            value.pop_back();
        }
        return *this;
    }

    // a, b >= 0, b > a; a-b = -(b-a), -a - (-b) = b-a = -(-b - (-a))
    else if (gt_abs(n, *this)) {
        bigint temp(n);
        temp -= *this;
        return *this = -temp;
    }

    // a == b
    else {
        value.clear();
        value.emplace_back(0);
        sign = POSITIVE;
        return *this;
    }
}

bigint& bigint::operator*=(const bigint& n) {
    if (!(*this) || !n) {
        return *this = bigint(0);
    }

    bigint p(0);
    p.value.resize(value.size() + n.value.size(), 0); // avoids for (...) if (...) 
    p.sign = static_cast<Sign>(sign * n.sign);

    for (int i = 0; i < n.value.size(); ++i) {
        uint64 carry = 0;
        for (int j = 0; j < value.size(); ++j) {
            uint64 new_val = carry + p.value[i+j] + value[j] * n.value[i];
            p.value[i+j] = new_val & mask; // new_val % base
            carry = new_val >> b_exp; // new_val / base
        }
        p.value[i+value.size()] = carry;
    }

    while (p.value.size() > 1 && p.value.back() == 0) {
        p.value.pop_back();
    }

    return *this = std::move(p);
}

// TODO
bigint& bigint::operator/=(const bigint& n) {
    return *this;
}

bigint& bigint::operator++() {
    return *this += 1;
}

bigint& bigint::operator++(int n) {
    return *this += 1;
}

bigint& bigint::operator--() {
    return *this -= 1;
}

bigint& bigint::operator--(int n) {
    return *this -= 1;
}

// TODO
bigint& bigint::operator%=(const bigint& n) {
    return *this;
}

bigint bigint::operator+(const bigint& n) const {
    return bigint(*this) += n;
}

bigint bigint::operator-(const bigint& n) const {
    return bigint(*this) -= n;
}

bigint bigint::operator*(const bigint& n) const {
    return bigint(*this) *= n;
}

bigint bigint::operator/(const bigint& n) const {
    return bigint(*this) /= n;
}

// TODO
bigint bigint::operator%(const bigint& n) const {
    return *this;
}

bigint bigint::operator-() const {
    bigint neg(*this);
    if (*this == 0) {
        neg.sign = POSITIVE;
        return neg;
    }
    neg.sign = sign == POSITIVE ? NEGATIVE : POSITIVE;
    return neg;
}

// returns abs(a) > abs(b)
bool bigint::gt_abs(const bigint& a, const bigint& b) {
    if (a.value.size() != b.value.size()) {
        return a.value.size() > b.value.size();
    }
    else {
        for (int i = a.value.size() - 1; i >= 0 ; --i) {
            if (a.value[i] != b.value[i]) {
                return a.value[i] > b.value[i];
            }
        }
    }
    return false;
}

bool bigint::operator>(const bigint& n) const {
    if (sign != n.sign) {
        return sign > n.sign;
    }
    else if (sign == POSITIVE) {
        return gt_abs(*this, n);
    }
    else {
        return gt_abs(n, *this);
    }
}

bool bigint::operator<(const bigint& n) const {
    return n > *this;
}

bool bigint::operator>=(const bigint& n) const {
    return !(n > *this);
}

bool bigint::operator<=(const bigint& n) const {
    return !(*this > n);
}

bool bigint::operator==(const bigint& n) const {
    return !(*this > n) && !(n > *this);
}

bool bigint::operator!=(const bigint& n) const {
    return !(*this == n);
}

bool bigint::operator!() const {
    return *this == 0;
}

// converts int/uint64/int64 n >= 0 -> bigint
void bigint::convert(uint64 n) {
    if (!n) {
        value.emplace_back(0);
    }
    while (n) {
        uint64 r = n & mask;
        n >>= b_exp;
        value.emplace_back(r);
    }
}

bigint bigint::abs() const {
    bigint abs(*this);
    abs.sign = POSITIVE;
    return abs;
}

// TODO base 10
std::string bigint::tostring(int str_len) const {
    if (str_len < 0) {
        return "Negative str_len\n";
    }
    static const char* hex_c = "0123456789abcdef";
    int shift = 4;
    std::string s{};
    s.reserve(value.size() * b_exp / shift);
    if (*this == bigint(0)) {
        return s + "0x0"; 
    }
    if (sign == NEGATIVE) {
        s += "-";
        if (str_len) {
            ++str_len;
        }
    }
    s += "0x";
    bool leading = true;
    for (int i = value.size() - 1; i >= 0; --i) {
        for (int j = b_exp - shift; j >= 0; j -= shift) {
            char c = hex_c[value[i] >> j & 0xf];
            if (leading && c == '0') {
                continue;
            }
            leading = false;
            s += c;
        }
    }
    if (str_len && s.length() >= str_len + 2) { 
        return s.substr(0, str_len + 2);
    }
    return s;
}
