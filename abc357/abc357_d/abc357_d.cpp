#include <iostream>
#define MOD 998244353
using namespace std;

long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}

long long mod_inv(long long a, long long mod) {
    return mod_pow(a, mod - 2, mod);
}

long long solve(long long n) {
    string n_str = to_string(n);
    long long length = n_str.size();
    long long x = mod_pow(10, length, MOD);
    long long x_n = mod_pow(x, n, MOD);
    
    long long numerator = (x_n - 1 + MOD) % MOD;
    long long denominator = (x - 1 + MOD) % MOD;
    long long denominator_inv = mod_inv(denominator, MOD);
    
    long long ans = n % MOD * numerator % MOD * denominator_inv % MOD;
    return ans;
}

int main() {
    long long n;
    cin >> n;
    cout << solve(n) << endl;
    return 0;
}
