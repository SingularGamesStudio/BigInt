#include <limits.h>

#include <algorithm>
#include <cassert>
#include <compare>
#include <iostream>
#include <string>
#include <vector>

using std::vector, std::string;

enum signs {
    neg = -1,
    zero = 0,
    pos = 1
};

const int MOD = 10;
class BigInteger {
   private:
    vector<int> data;
    signs sign;

    int get(size_t id) const {
        if (id < data.size())
            return data[id];
        return 0;
    }

    void upgrade(int x, int& zeros) {
        if (x) {
            for (int j = 0; j < zeros; j++) {
                data.push_back(0);
            }
            data.push_back(x);
            zeros = 0;
        } else
            ++zeros;
    }

    BigInteger abs() const {
        BigInteger res = BigInteger(*this);
        if (res.sign == signs::neg)
            res.sign = signs::pos;
        return res;
    }

    void add1() {
        int delta = 1;
        for (size_t i = 0; i < data.size(); i++) {
            data[i] += delta;
            delta = data[i] / MOD;
            data[i] = data[i] % MOD;
            if (delta == 0)
                break;
        }
        if (delta != 0)
            data.push_back(delta);
    }

    void sub1() {
        int delta = -1;
        for (size_t i = 0; i < data.size(); i++) {
            data[i] += delta;
            if (data[i] < 0) {
                delta = -1;
                data[i] += MOD;
            } else
                break;
        }
        while (data.size() && data[data.size() - 1] == 0)
            data.pop_back();
        data.shrink_to_fit();
    }

    BigInteger add(const BigInteger& second, signs sign1 = signs::pos) const {
        if (second.sign == signs::zero || sign1 == signs::zero)
            return BigInteger(*this);
        if (sign == signs::zero)
            return BigInteger(second);
        BigInteger res = BigInteger();
        int size = std::max(data.size(), second.data.size()) + 1;
        res.data.reserve(size);
        int delta = 0;
        int zeros = 0;
        signs sign2 = ((sign1 == signs::neg && second.sign == signs::pos) || (sign1 == signs::pos && second.sign == signs::neg)) ? signs::neg : signs::pos;
        if (sign == sign2) {
            for (int i = 0; i < size; i++) {
                int now = get(i) + second.get(i) + delta;
                delta = now / MOD;
                now = now % MOD;
                res.upgrade(now, zeros);
            }
            res.sign = sign;
        } else {
            BigInteger abs1 = abs();
            BigInteger abs2 = second.abs();
            bool swapped = false;
            if (abs1 < abs2) {
                std::swap(abs1, abs2);
                swapped = true;
            }
            for (int i = 0; i < size; i++) {
                int now = abs1.get(i) - abs2.get(i) + delta;
                if (now < 0) {
                    delta = -1;
                    now = now + MOD;
                } else
                    delta = 0;
                res.upgrade(now, zeros);
            }
            if (zeros == size) {
                res.sign = zero;
            } else if (!swapped)
                res.sign = sign;
            else
                res.sign = sign2;
        }
        res.data.shrink_to_fit();
        return res;
    }

    BigInteger mul(const BigInteger& second) const {
    }

   public:
    BigInteger() : sign(signs::zero) {}

    BigInteger(long long n) {
        if (n == 0) {
            sign = signs::zero;
        } else if (n > 0) {
            sign = signs::pos;
        } else
            sign = signs::neg;
        n = std::abs(n);
        while (n > 0) {
            data.push_back(n % MOD);
            n /= MOD;
        }
    }

    BigInteger(unsigned long long n) {
        if (n == 0) {
            sign = signs::zero;
        } else
            sign = signs::pos;
        while (n > 0) {
            data.push_back(n % MOD);
            n /= MOD;
        }
    }

    BigInteger(const BigInteger& base) {
        data = base.data;
        sign = base.sign;
    }

    std::strong_ordering operator<=>(const BigInteger& second) const {
        if (sign == signs::pos && second.sign != signs::pos)
            return std::strong_ordering::greater;
        if (sign == signs::neg && second.sign != signs::neg)
            return std::strong_ordering::less;
        if (sign == signs::zero)
            return 0 <=> static_cast<int>(second.sign);
        const BigInteger& swapped1 = (sign == signs::pos) ? *this : second;
        const BigInteger& swapped2 = (sign == signs::pos) ? second : *this;
        if (swapped1.data.size() != swapped2.data.size())
            return swapped1.data.size() <=> swapped2.data.size();
        for (int i = data.size() - 1; i >= 0; --i) {
            if (data[i] != second.data[i])
                return swapped1.data[i] <=> swapped2.data[i];
        }
        return std::strong_ordering::equal;
    }

    BigInteger& operator++() {
        if (sign == signs::zero) {
            sign = signs::pos;
            data.push_back(1);
        } else if (sign == signs::pos) {
            add1();
        } else {
            sub1();
            if (data.size() == 0)
                sign = signs::zero;
        }
        return *this;
    }

    BigInteger& operator--() {
        if (sign == signs::zero) {
            sign = signs::neg;
            data.push_back(1);
        } else if (sign == signs::neg) {
            add1();
        } else {
            sub1();
            if (data.size() == 0)
                sign = signs::zero;
        }
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger res = BigInteger(*this);
        ++(*this);
        return res;
    }

    BigInteger operator--(int) {
        BigInteger res = BigInteger(*this);
        --(*this);
        return res;
    }

    BigInteger operator-() const {
        BigInteger res = BigInteger(*this);
        if (sign == signs::pos)
            res.sign = signs::neg;
        else if (sign == signs::neg)
            res.sign = signs::pos;
        return res;
    }

    string toString() const {
        string s = "";
        if (sign == signs::neg)
            s += '-';
        for (int i = data.size() - 1; i >= 0; i--) {
            s += '0' + data[i];
        }
        if (sign == signs::zero)
            s = "0";
        return s;
    }

    explicit operator bool() const {
        return sign != signs::zero;
    }

    friend BigInteger operator+(const BigInteger& first, const BigInteger& second);
    friend BigInteger operator-(const BigInteger& first, const BigInteger& second);
    friend std::ostream& operator<<(std::ostream& output, const BigInteger& val);
    friend std::istream& operator>>(std::istream& input, BigInteger& val);
};

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    return first.add(second);
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    return first.add(second, signs::neg);
}

std::ostream& operator<<(std::ostream& output, const BigInteger& val) {
    output << val.toString();
    return output;
}

std::istream& operator>>(std::istream& input, BigInteger& val) {
    val.data.clear();
    val.sign = signs::pos;
    char c;
    while (input.peek()) {
        if (std::isspace(input.peek()))
            input.get(c);
        else
            break;
    }
    bool started = false;
    while (input.get(c) && !std::isspace(c)) {
        if (started) {
            assert(c >= '0' && c <= '9');
            val.data.push_back(c - '0');
        } else {
            if (c == '-')
                val.sign = signs::neg;
            if (c > '0' && c <= '9') {
                started = true;
                val.data.push_back(c - '0');
            }
        }
    }
    if (!started)
        val.sign = signs::zero;
    return input;
}

BigInteger operator"" _bi(unsigned long long x) {
    BigInteger res(x);
    return x;
}

/*
    умножение, деление, остаток по модулю, составное присваивание с этими операциями.
*/