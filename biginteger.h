#include <limits.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <compare>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

using std::vector, std::string, std::deque;

enum signs {
    neg = -1,
    zero = 0,
    pos = 1
};

signs mulsigns(signs a, signs b) {
    if (a == signs::zero || b == signs::zero)
        return signs::zero;
    if (a == b)
        return signs::pos;
    return signs::neg;
}

int reversebits(int x, int pw2) {
    int res = 0;
    for (int i = 0; i < pw2; i++) {
        if ((x >> i) & 1)
            res |= (1 << (pw2 - i - 1));
    }
    return res;
}

const int MOD = 1000;
const int DIGITS = 3;

class Complex {
   public:
    double r;
    double i;

    Complex() : r(0), i(0) {}

    Complex(double x) : r(x), i(0) {}

    Complex(double r, double i) : r(r), i(i) {}

    Complex(const Complex &tocopy) : r(tocopy.r), i(tocopy.i) {}

    Complex operator+(const Complex &second) const {
        return Complex(r + second.r, i + second.i);
    }

    Complex operator-(const Complex &second) const {
        return Complex(r - second.r, i - second.i);
    }

    Complex operator*(const Complex &second) const {
        return Complex(r * second.r - i * second.i, i * second.r + r * second.i);
    }

    Complex &operator*=(const Complex &second) {
        double temp = r;
        double temp1 = second.r;
        r = r * second.r - i * second.i;
        i = i * temp1 + temp * second.i;
        return *this;
    }

    Complex &operator=(const Complex &second) {
        r = second.r;
        i = second.i;
        return *this;
    }
};

std::ostream &operator<<(std::ostream &output, const Complex &val) {
    output << ((double)(int)(val.r * 100.0)) / 100 << "+" << ((double)(int)(val.i * 100.0)) / 100 << "i";
    return output;
}

void FFT(Complex *a, int n, Complex q) {
    int pw2 = 0;
    while ((1 << pw2) < n)
        pw2++;
    for (int i = 0; i < n; i++) {
        int rev = reversebits(i, pw2);
        if (i < rev)
            std::swap(a[i], a[rev]);
    }
    for (int l = 2; l <= n; l *= 2) {
        Complex cur = q;
        for (int ll = n; ll > l; ll /= 2)
            cur *= cur;
        for (int start = 0; start < n; start += l) {
            int mid = start + l / 2;
            Complex qdeg = 1;
            int pos = start;
            while (pos < mid) {
                Complex u = Complex(a[pos]);
                Complex v = Complex(a[pos + l / 2]) * qdeg;
                a[pos] = u + v;
                a[pos + l / 2] = u - v;
                qdeg *= cur;
                pos++;
            }
        }
    }
}

class PoweredInteger;

class BigInteger {
   private:
    int get(size_t id) const {
        if (id < data.size()) return data[id];
        return 0;
    }

    void upgrade(int x, int &zeros) {
        if (x) {
            for (int j = 0; j < zeros; j++) {
                data.push_back(0);
            }
            data.push_back(x);
            zeros = 0;
        } else
            ++zeros;
    }

    void add1() {  // add 1 to the abs
        if (sign == signs::zero) {
            data.push_back(1);
            sign = signs::pos;
            return;
        }
        int delta = 1;
        for (size_t i = 0; i < data.size(); i++) {
            data[i] += delta;
            delta = data[i] / MOD;
            data[i] = data[i] % MOD;
            if (delta == 0) break;
        }
        if (delta != 0) data.push_back(delta);
    }

    void sub1() {  // substract 1 from the abs
        if (sign == signs::zero) {
            data.push_back(1);
            sign = signs::neg;
            return;
        }
        int delta = -1;
        for (size_t i = 0; i < data.size(); i++) {
            data[i] += delta;
            if (data[i] < 0) {
                delta = -1;
                data[i] += MOD;
            } else
                break;
        }
        while (data.size() && data[data.size() - 1] == 0) data.pop_back();
        data.shrink_to_fit();
        if (data.size() == 0) sign = signs::zero;
    }

    void swap(BigInteger &other) {
        std::swap(data, other.data);
        std::swap(sign, other.sign);
    }

    std::strong_ordering cmpabs(const BigInteger &second) const {
        if (data.size() != second.data.size())
            return data.size() <=> second.data.size();
        for (int i = data.size() - 1; i >= 0; --i) {
            if (data[i] != second.data[i])
                return data[i] <=> second.data[i];
        }
        return std::strong_ordering::equal;
    }

    BigInteger add(const BigInteger &second, signs sign1 = signs::pos) const {
        if (second.sign == signs::zero || sign1 == signs::zero)
            return BigInteger(*this);
        if (sign == signs::zero) return BigInteger(second);
        BigInteger res = BigInteger();
        int size = std::max(data.size(), second.data.size()) + 1;
        int delta = 0;
        int zeros = 0;
        signs sign2 = mulsigns(second.sign, sign1);
        if (sign == sign2) {  // sum of modules
            for (int i = 0; i < size; i++) {
                int now = get(i) + second.get(i) + delta;
                delta = now / MOD;
                now = now % MOD;
                res.upgrade(now, zeros);
            }
            res.sign = sign;
        } else {  // diff of modules
            const BigInteger *abs1 = this;
            const BigInteger *abs2 = &second;
            bool swapped = false;
            if (cmpabs(second) == std::strong_ordering::less) {
                std::swap(abs1, abs2);
                swapped = true;
            }
            for (int i = 0; i < size; i++) {
                int now = (*abs1).get(i) - (*abs2).get(i) + delta;
                if (now < 0) {
                    delta = -1;
                    now = now + MOD;
                } else
                    delta = 0;
                res.upgrade(now, zeros);
            }
            if (zeros == size) {
                res.sign = zero;
            } else
                res.sign = swapped ? sign2 : sign;
        }
        return res;
    }

    BigInteger mul(const BigInteger &second) const {
        if (sign == signs::zero || second.sign == signs::zero)
            return 0;
        size_t n = 1;
        while (n < data.size() || n < second.data.size())
            n *= 2;
        n *= 2;
        Complex *a = new Complex[n]();
        Complex *b = new Complex[n]();
        for (size_t i = 0; i < data.size(); i++)
            a[i] = data[i];
        for (size_t i = 0; i < second.data.size(); i++)
            b[i] = second.data[i];
        double phi = 2 * acos(-1) / static_cast<double>(n);
        Complex q = Complex(cos(phi), sin(phi));
        FFT(a, n, q);
        FFT(b, n, q);
        for (size_t i = 0; i < n; i++)
            a[i] *= b[i];
        FFT(a, n, Complex(cos(-phi), sin(-phi)));
        BigInteger res;
        long long delta = 0;
        size_t pos = 0;
        int cntzero = 0;
        while (pos < n || delta) {
            delta += round(a[pos].r / n);
            if (delta % MOD) {
                for (int i = 0; i < cntzero; i++)
                    res.data.push_back(0);
                cntzero = 0;
                res.data.push_back(delta % MOD);
            } else
                cntzero++;
            delta /= MOD;
            pos++;
        }
        res.sign = mulsigns(sign, second.sign);
        delete[] a;
        delete[] b;
        return res;
    }

    long double getfirst(size_t cnt) const {
        long double res = 0;
        long double power = static_cast<long double>(1) / MOD;
        for (size_t i = 0; i < std::min(cnt / DIGITS + 2, data.size()); i++) {
            res += power * data[data.size() - i - 1];
            power /= MOD;
        }
        while (res < 0.1)
            res *= 10;
        return res;
    }

   public:
    deque<int> data;
    signs sign;

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

    BigInteger(const BigInteger &base) {
        data = base.data;
        sign = base.sign;
    }

    BigInteger &operator=(const BigInteger &second) {
        data = second.data;
        sign = second.sign;
        return *this;
    }

    BigInteger &operator++() {
        if (sign == signs::zero) {
            sign = signs::pos;
            data.push_back(1);
        } else if (sign == signs::pos) {
            add1();
        } else {
            sub1();
        }
        return *this;
    }

    BigInteger &operator--() {
        if (sign == signs::zero) {
            sign = signs::neg;
            data.push_back(1);
        } else if (sign == signs::neg) {
            add1();
        } else {
            sub1();
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
        res.sign = mulsigns(res.sign, signs::neg);
        return res;
    }

    string toString() const {
        string s = "";
        if (sign == signs::neg) s += '-';
        for (int i = data.size() - 1; i >= 0; i--) {
            string s0 = std::to_string(data[i]);
            if (i != static_cast<int>(data.size()) - 1)
                while (s0.size() < 3) {
                    s0 = "0" + s0;
                }
            s += s0;
        }
        if (sign == signs::zero) s = "0";
        return s;
    }

    explicit operator bool() const { return sign != signs::zero; }

    BigInteger &operator+=(const BigInteger &second);
    BigInteger &operator-=(const BigInteger &second);
    BigInteger &operator*=(const BigInteger &second);
    BigInteger &operator/=(const BigInteger &second);
    BigInteger &operator%=(const BigInteger &second);

    friend BigInteger operator+(const BigInteger &first, const BigInteger &second);
    friend BigInteger operator-(const BigInteger &first, const BigInteger &second);
    friend std::ostream &operator<<(std::ostream &output, const BigInteger &val);
    friend std::istream &operator>>(std::istream &input, BigInteger &val);
    friend BigInteger operator*(const BigInteger &first, const BigInteger &second);
    friend BigInteger operator/(const BigInteger &first, const BigInteger &second);
    friend std::strong_ordering operator<=>(const BigInteger &first, const BigInteger &second);
    friend bool operator==(const BigInteger &first, const BigInteger &second);
    friend PoweredInteger divide(const BigInteger &first, const BigInteger &second, size_t precision);
    friend BigInteger gcd(BigInteger a, BigInteger b);
};

std::strong_ordering operator<=>(const BigInteger &first, const BigInteger &second) {
    if (first.sign == signs::pos && second.sign != signs::pos)
        return std::strong_ordering::greater;
    if (first.sign == signs::neg && second.sign != signs::neg)
        return std::strong_ordering::less;
    if (first.sign == signs::zero) return 0 <=> static_cast<int>(second.sign);
    const BigInteger &swapped1 = (first.sign == signs::pos) ? first : second;
    const BigInteger &swapped2 = (first.sign == signs::pos) ? second : first;
    return swapped1.cmpabs(swapped2);
}

bool operator==(const BigInteger &first, const BigInteger &second) {
    if (first.sign != second.sign)
        return false;
    if (first.data.size() != second.data.size())
        return false;
    for (int i = first.data.size() - 1; i >= 0; --i) {
        if (first.data[i] != second.data[i])
            return false;
    }
    return true;
}

bool operator!=(const BigInteger &first, const BigInteger &second) {
    return !(first == second);
}

BigInteger operator+(const BigInteger &first, const BigInteger &second) {
    return first.add(second);
}

BigInteger operator*(const BigInteger &first, const BigInteger &second) {
    return first.mul(second);
}

BigInteger operator-(const BigInteger &first, const BigInteger &second) {
    return first.add(second, signs::neg);
}

class PoweredInteger {
   public:
    BigInteger val;
    long long power;
    PoweredInteger(BigInteger val, long long power) : val(val), power(power) {}
    PoweredInteger() : val(BigInteger(0)), power(0) {}
    PoweredInteger(const PoweredInteger &tocopy) : val(tocopy.val), power(tocopy.power) {}

    PoweredInteger &operator=(const PoweredInteger &second) {
        val = second.val;
        power = second.power;
        return *this;
    }

    PoweredInteger operator*(const PoweredInteger &second) const {
        PoweredInteger res;
        res.val = val * second.val;
        res.power = power + second.power;
        return res;
    }

    PoweredInteger operator+(const PoweredInteger &second) const {
        long long minpower = std::min(power, second.power);
        PoweredInteger a = PoweredInteger(*this);
        PoweredInteger b = PoweredInteger(second);
        if (a.val.sign != signs::zero)
            while (a.power > minpower) {
                a.val.data.push_front(0);
                a.power--;
            }
        if (b.val.sign != signs::zero)
            while (b.power > minpower) {
                b.val.data.push_front(0);
                b.power--;
            }
        return PoweredInteger(a.val + b.val, minpower);
    }

    PoweredInteger operator-(const PoweredInteger &second) const {
        long long minpower = std::min(power, second.power);
        PoweredInteger a = PoweredInteger(*this);
        PoweredInteger b = PoweredInteger(second);
        if (a.val.sign != signs::zero)
            while (a.power > minpower) {
                a.val.data.push_front(0);
                a.power--;
            }
        if (b.val.sign != signs::zero)
            while (b.power > minpower) {
                b.val.data.push_front(0);
                b.power--;
            }
        return PoweredInteger(a.val - b.val, minpower);
    }

    void cut(int length) {  // rounds to closest
        int delta = 0;
        while (val.data.size() > static_cast<size_t>(length)) {
            delta = val.data.front();
            val.data.pop_front();
            power++;
        }
        if (delta >= MOD / 2)
            val++;
    }

    string toString(int precision10) const {  // TODO: rewrite this shit
        string s = val.toString();
        string o = "";
        for (int j = 0; j < DIGITS; j++)
            o += '0';
        if (power >= 0) {
            for (int i = 0; i < power; i++)
                s += o;
            if (precision10 > 0) {
                s.push_back('.');
                for (int i = 0; i < precision10; i++)
                    s.push_back('0');
            }
            return s;
        }

        bool neg = (s[0] == '-');
        string s0 = "";
        for (size_t i = neg; i < s.size(); i++)
            s0 += s[i];
        while (static_cast<int>(s0.length()) < -DIGITS * power + 1) {
            s0 = "0" + s0;
        }
        string s1 = neg ? "-" : "";
        int cnt = 0;
        for (size_t i = 0; i < s0.length(); i++) {
            s1.push_back(s0[i]);
            if (cnt) {
                cnt++;
                if (cnt >= precision10 + 1)
                    break;
            }
            if (static_cast<int>(s0.length()) - static_cast<int>(i) - 1 == -DIGITS * power) {
                if (precision10 == 0)
                    break;
                s1.push_back('.');
                cnt++;
            }
        }
        while (cnt < precision10 + 1) {
            s1.push_back(0);
            cnt++;
        }
        return s1;
    }

    explicit operator BigInteger() {  // rounds to zero
        BigInteger res = val;
        if (power > 0) {
            for (int i = 0; i < power; i++) {
                res.data.push_front(0);
            }
        } else {
            for (int i = 0; i < -power; i++) {
                if (res.data.size() == 0)
                    return BigInteger(0);
                res.data.pop_front();
            }
            if (res.data.size() == 0) {
                return BigInteger(0);
            }
        }
        return res;
    }

    explicit operator double() {
        double pw10 = pow(MOD, power);
        double ans = 0;
        for (size_t i = 0; i < val.data.size(); i++) {
            ans += pw10 * val.data[i];
            pw10 *= MOD;
        }
        return ans;
    }
};

std::ostream &operator<<(std::ostream &output, const BigInteger &val) {
    output << val.toString();
    return output;
}

std::ostream &operator<<(std::ostream &output, const PoweredInteger &val) {
    output << val.toString(100);
    return output;
}

PoweredInteger divide(const BigInteger &first, const BigInteger &second, size_t precision) {
    assert(second.sign != signs::zero);
    if (first.sign == signs::zero)
        return PoweredInteger(BigInteger(0), 0);
    signs ressign = mulsigns(first.sign, second.sign);
    long double firstiter = 1.0 / second.getfirst(5);
    PoweredInteger rev = PoweredInteger(BigInteger(static_cast<long long>(1000000000.0 * firstiter)), -9 / DIGITS);
    PoweredInteger two = PoweredInteger(2, 0);
    int sizediff = 1;
    int temp = second.data[second.data.size() - 1];
    while (temp * 10 < MOD) {
        temp *= 10;
        sizediff *= 10;
    }
    PoweredInteger cur = PoweredInteger(second * sizediff, -second.data.size());
    cur.val.sign = signs::pos;
    rev.val.sign = signs::pos;
    int maxsigns = 1;
    while (static_cast<size_t>(maxsigns) < first.data.size() + precision) {
        maxsigns *= 2;
    }
    maxsigns *= 4;
    for (int iter = 1; static_cast<size_t>((1 << std::max(0, (iter - 3)))) < first.data.size() + precision; iter++) {
        // std::cout << rev << "\n";
        rev = rev * (two - cur * rev);
        rev.cut(maxsigns);
    }
    // std::cout << rev << "\n";
    // std::cout << cur << "\n";
    cur = PoweredInteger(first * sizediff, -second.data.size());
    // std::cout << cur << "\n";
    rev = cur * rev;
    // std::cout << rev << "\n";
    rev.val.sign = ressign;
    return rev;
}

BigInteger operator/(const BigInteger &first, const BigInteger &second) {
    PoweredInteger res = divide(first, second, 0);
    signs ressign = res.val.sign;
    BigInteger ans = static_cast<BigInteger>(res);
    res.val.data.clear();
    ans.add1();
    if ((ans * second).cmpabs(first) == std::strong_ordering::greater)  // for example, 228/3 = 75.999999999999999999999999999999=76, but rounded down to 75
        ans.sub1();
    if (ans.sign != signs::zero)
        ans.sign = ressign;
    return ans;
}

BigInteger operator%(const BigInteger &first, const BigInteger &second) {
    return first - (first / second) * second;
}

std::istream &operator>>(std::istream &input, BigInteger &val) {
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
    vector<int> pool;
    while (input.get(c) && !std::isspace(c)) {
        if (started) {
            assert(c >= '0' && c <= '9');
            pool.push_back(c - '0');
        } else {
            if (c == '-') val.sign = signs::neg;
            if (c > '0' && c <= '9') {
                started = true;
                pool.push_back(c - '0');
            }
        }
    }
    std::reverse(pool.begin(), pool.end());
    for (int i = 0; i < static_cast<int>(pool.size()); i += DIGITS) {
        int x = 0;
        for (int j = i + DIGITS - 1; j >= i; j--) {
            if (j < static_cast<int>(pool.size())) {
                x *= 10;
                x += pool[j];
            }
        }
        val.data.push_back(x);
    }
    if (!started) val.sign = signs::zero;
    return input;
}

BigInteger operator"" _bi(unsigned long long x) {
    BigInteger res(x);
    return x;
}

BigInteger &BigInteger::operator+=(const BigInteger &second) {
    *this = *this + second;
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &second) {
    *this = *this - second;
    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &second) {
    *this = *this * second;
    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &second) {
    *this = *this / second;
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &second) {
    *this = *this % second;
    return *this;
}

BigInteger gcd(BigInteger a, BigInteger b) {
    if (a.sign == signs::zero)
        return b;
    if (b.sign == signs::zero)
        return a;
    a.sign = signs::pos;
    b.sign = signs::pos;
    while (b.sign != signs::zero) {
        a.swap(b);
        b = b % a;
    }
    return a;
}

class Rational {
   public:
    BigInteger p;
    BigInteger q;
    Rational() : p(0), q(1) {}
    Rational(const BigInteger &x) : p(x), q(1) {}
    Rational(const BigInteger &first, const BigInteger &second) : p(first), q(second) {
        assert(q.sign != signs::zero);
        if (q.sign == signs::neg) {
            q.sign = signs::pos;
            p.sign = mulsigns(p.sign, signs::neg);
        }
    }
    Rational(long long x) : p(x), q(1) {}
    Rational(const Rational &tocopy) : p(tocopy.p), q(tocopy.q) {}

    Rational &operator=(const Rational &second) {
        p = second.p;
        q = second.q;
        return *this;
    }

    Rational operator-() const {
        Rational res = *this;
        res.p.sign = mulsigns(res.p.sign, signs::neg);
        return res;
    }

    void normalize() {
        BigInteger div = gcd(p, q);
        p /= div;
        q /= div;
    }

    string toString() {
        normalize();
        if (p.sign == signs::zero)
            return "0";
        if (q == BigInteger(1))
            return p.toString();
        return p.toString() + "/" + q.toString();
    }

    string asDecimal(size_t precision10 = 0) {
        PoweredInteger res = divide(p, q, precision10 / DIGITS + 2);
        return res.toString(precision10);
    }

    explicit operator double() {
        PoweredInteger res = divide(p, q, 20);
        return static_cast<double>(res);
    }

    Rational &operator+=(const Rational &second);
    Rational &operator-=(const Rational &second);
    Rational &operator*=(const Rational &second);
    Rational &operator/=(const Rational &second);
};

Rational operator+(const Rational &first, const Rational &second) {
    return Rational(first.p * second.q + second.p * first.q, first.q * second.q);
}

Rational operator-(const Rational &first, const Rational &second) {
    return Rational(first.p * second.q - second.p * first.q, first.q * second.q);
}

Rational operator*(const Rational &first, const Rational &second) {
    return Rational(first.p * second.p, first.q * second.q);
}

Rational operator/(const Rational &first, const Rational &second) {
    assert(second.p.sign != signs::zero);
    return Rational(first.p * second.q, first.q * second.p);
}

Rational &Rational::operator+=(const Rational &second) {
    *this = *this + second;
    return *this;
}

Rational &Rational::operator-=(const Rational &second) {
    *this = *this - second;
    return *this;
}

Rational &Rational::operator*=(const Rational &second) {
    *this = *this * second;
    return *this;
}

Rational &Rational::operator/=(const Rational &second) {
    *this = *this / second;
    return *this;
}

std::strong_ordering operator<=>(const Rational &first, const Rational &second) {
    if (first.p.sign == signs::pos && second.p.sign != signs::pos)
        return std::strong_ordering::greater;
    if (first.p.sign == signs::neg && second.p.sign != signs::neg)
        return std::strong_ordering::less;
    if (first.p.sign == signs::zero) return 0 <=> static_cast<int>(second.p.sign);
    return first.p * second.q <=> second.p * first.q;
}

bool operator==(const Rational &first, const Rational &second) {
    return first.p * second.q == second.p * first.q;
}

bool operator!=(const Rational &first, const Rational &second) {
    return !(first == second);
}