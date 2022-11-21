#include "biginteger.h"

#include <iostream>

using namespace std;

int main() {
    BigInteger a;
    BigInteger b;
    cin >> a >> b;
    // cout << a << " " << b << " " << BigInteger(455) << endl;
    cout << (a * b) << "\n";
}