#include "biginteger.h"

#include <bits/stdc++.h>

using namespace std;

bool st10(long long x) {
    while (x % 10 == 0)
        x /= 10;
    return x == 1;
}

int main() {
    BigInteger a;
    BigInteger b;
    mt19937_64 rnd(42);
    // cin >> a >> b;
    // cout << a / b;
    int cnt = 0;
    while (1) {
        cnt++;
        if (cnt >= 20) {
            cnt = 0;
            cout << "20 epochs done" << endl;
        }
        long long x = rnd() % (long long)1e17, y = rnd() % (long long)1e17;
        if (!st10(x) && !st10(y)) {
            a = BigInteger(x), b = BigInteger(y);
            if (x / y != (a / b)) {
                cout << x << " " << y << endl;
                cout << "true: " << setprecision(20) << x / y << "\nfalse: " << a / b << endl;
                return 0;
            }
        }
    }
    return 0;
}