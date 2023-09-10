#include <bits/stdc++.h>
#define OVERLOAD_REP(_1, _2, _3, name, ...) name
#define REP1(i, n) for (auto i = std::decay_t<decltype(n)>{}; (i) != (n); ++(i))
#define REP2(i, l, r) for (auto i = (l); (i) != (r); ++(i))
#define rep(...) OVERLOAD_REP(__VA_ARGS__, REP2, REP1)(__VA_ARGS__)
#define OVERLOAD_RREP(_1, _2, _3, name, ...) name
#define RREP1(i, n) for (auto i = (n); (i) != (std::decay_t<decltype(n)>{}); --(i))
#define RREP2(i, l, r) for (auto i = (l); (i) != (r); --(i))
#define rrep(...) OVERLOAD_RREP(__VA_ARGS__, RREP2, RREP1)(__VA_ARGS__)
#define all(x) (x).begin,(x).end()
#define rall(x) (x).rbegin, (x).rend
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
#define YESNO(T) if(T){cout<<"YES"<<endl;}else{cout<<"NO"<<endl;}
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
using namespace std;
using ll=long long;
using ld=long double;
using ull=unsigned long long;
using vi = vector<int>;
using vl = vector<long>;
using vll = vector<long long>;
using vvi = vector<vi>;
using vvl = vector<vl>;
using pii = pair<int, int>;
using pli = pair<ll, int>;
using pll = pair<ll, ll>;
using Graph = vector<vector<ll>>;
const string ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const string abc="abcdefghijklmnopqrstuvwxyz";
#define gcd __gcd
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define print(n) cout << n << endl
template <typename T>
inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false)); }
template <typename T>
inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false)); }
const ll mod=998244353;
const ll MOD=1000000007;
const ld PI=3.1415265358979323;
//数学系
//繰り返し二乗法
ll modpow(ll a, ll b, ll mod) {
    if(a==0) {
        return 0;
    }
    ll res = 1;
    while(b > 0) {
        if(b & 1) res = res * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return res;
}
    // 拡張 Euclid の互除法
    long long extGCD(long long a, long long b, long long& p, long long& q) {
        if (b == 0) { p = 1; q = 0; return a; }
        long long d = extGCD(b, a % b, q, p);
        q -= a / b * p;
        return d;
    }
    long long extGcd(long long a, long long b, long long& p, long long& q) {
        if (b == 0) { p = 1; q = 0; return a; }
        long long d = extGcd(b, a % b, q, p);
        q -= a / b * p;
        return d;
    }
//逆元
ll modinv(ll a, ll m) {
    ll b = m, u = 1, v = 0;
    while (b) {
        ll t = a / b;
        a -= t * b; swap(a, b);
        u -= t * v; swap(u, v);
    }
    u %= m; 
    if (u < 0) u += m;
    return u;
}
//割り算
ll div(ll a, ll b, ll mod) {
return a * modinv(b,mod) % mod;
}
//階乗
vll fact(1);
void fac(ll n, ll m) {
    fact[0]=1;
    for(int i=1;i<=n;i++) {
        fact.pb(fact[i-1]*i);
        fact[i]%=m;
    }
}
//組み合わせ
ll comb(ll n,ll k,ll mod) {
return div(fact[n],(fact[k]*fact[n-k] % mod),mod);
}

int main() {
    int n;
    cin>>n;
    int ans=0;
    for(int i=1;i<=n;i+=2) {
        int a=0;
        for(int j=1;j<=i;j++){
            if(i%j==0) a++;
        }
        if(a==8) ans++;
    }
    print(ans);
}