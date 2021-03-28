#pragma once 

#include <string>
#include <vector>
#include <istream>
#include <ostream>

class bigint {
    private:
        using uint32 = std::uint32_t;
        using int64 = std::int64_t;
        using uint64 = std::uint64_t;

        enum Sign {
            NEGATIVE = -1,
            POSITIVE = 1
        };

        static constexpr uint64 b_exp = 32;
        static constexpr uint64 base = 0x1ULL << b_exp;
        static constexpr uint64 mask = base - 1;

        Sign sign;
        std::vector<uint64> value{};

        void convert(uint64 n);
        static bool gt_abs(const bigint& a, const bigint& b);

    public:
        bigint(int n = 0);
        bigint(int64 n);
        bigint(uint64 n);
        bigint(const bigint& n);
        bigint(std::string s);

        bigint& operator=(int n);
        bigint& operator=(int64 n);
        bigint& operator=(const bigint& n);
        bigint& operator=(bigint&& n);

        bigint operator-() const;

        bigint& operator+=(const bigint& n);
        bigint& operator-=(const bigint& n);
        bigint& operator*=(const bigint& n);
        bigint& operator/=(const bigint& n);
        bigint& operator%=(const bigint& n);

        bigint& operator++();
        bigint& operator++(int n);
        bigint& operator--();
        bigint& operator--(int n);

        bigint operator+(const bigint& n) const;
        bigint operator-(const bigint& n) const;
        bigint operator*(const bigint& n) const;
        bigint operator/(const bigint& n) const;
        bigint operator%(const bigint& n) const;

        bool operator>(const bigint& n) const;
        bool operator<(const bigint& n) const;
        bool operator>=(const bigint& n) const;
        bool operator<=(const bigint& n) const;
        bool operator==(const bigint& n) const;
        bool operator!=(const bigint& n) const;
        bool operator!() const;

        bigint abs() const;
        std::string tostring(int str_len = 0) const;
        friend std::istream& operator>>(std::istream& in, bigint& n);
        friend std::ostream& operator<<(std::ostream out, const bigint& n);
};
