#include <bits/stdc++.h>

using namespace std;

int main() {
    int n,m,p;
    cin >> n >> m >> p;
    if (n>=m) cout << (n-m) / p + 1 << endl;
    else cout << 0 << endl;
}