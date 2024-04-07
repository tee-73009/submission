#line 1 "code.cpp"
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#include <bits/stdc++.h>
using namespace std;
#define OVERLOAD_3(_1, _2, _3, name, ...) name
#define REP1(i, n) for (auto i = std::decay_t<decltype(n)>{}; (i) != (n); ++(i))
#define REP2(i, l, r) for (auto i = (l); (i) != (r); ++(i))
#define rep(...) OVERLOAD_3(__VA_ARGS__, REP2, REP1)(__VA_ARGS__)
signed main() {
    ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    long long n,m;
    cin >> n >> m;
    vector<vector<long long>> a(n,vector<long long>(m));
    rep(i,n) rep(j,m) cin >> a[i][j];
    long long ans=0;
    rep(i,n){
        rep(j,i+1,n){
            long long b=0;
            rep(k,m){
                b^=(a[i][k]==a[j][k]);
            }
            ans+=b;
        }
    }
    cout << ans << endl;
}
