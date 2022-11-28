#include "biginteger.h"

#include <bits/stdc++.h>

using namespace std;

int main() {
    BigInteger a;
    BigInteger b;
    mt19937_64 rnd(42);
    // int cnt = 0;
    cout << setprecision(30) << (Rational(1) / Rational(-25)).asDecimal(4);
    /*while (1) {
        cnt++;
        if (cnt >= 20) {
            cnt = 0;
            cout << "20 epochs done" << endl;
        }
        long long x = rnd() % (long long)1e9, y = rnd() % (long long)1e9;
        a = BigInteger(x), b = BigInteger(y);
        if ((double)x / y != (double)Rational(a, b)) {
            cout << x << " " << y << endl;
            cout << "true: " << setprecision(20) << (double)x / y << "\nfalse: " << (double)Rational(a, b) << endl;
            return 0;
        }
    }*/
    return 0;
}