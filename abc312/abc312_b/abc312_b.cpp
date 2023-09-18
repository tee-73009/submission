#ifdef INCLUDED_MAIN

int main() {
  int n,m;
  cin >> n >> m;
  vector<string> a(110);
  rep(i,n)cin >> a[i];
  for(int i=0;i+8<=n;i++){
    for(int j=0;j+8<=m;j++) {
      bool can=true;
      rep(ii,4){
        rep(jj,4){
          if(ii!=3&&jj!=3){
            if(a[i+ii][j+jj]!='#')can=false;
            if(a[i+8-ii][j+8-jj]!='#')can=false;
          }
          else{
            if(a[i+ii][j+jj]!='.')can=false;
            if(a[i+8-ii][j+8-jj]!='.')can=false;
          }
        }
      }
      if(can)cout << i+1 << ' ' << j+1 << endl;
    }
  }
}

#else
//ここからテンプレート
#include <bits/stdc++.h>
#define OVERLOAD_REP(_1, _2, _3, name, ...) name
#define REP1(i, n) for (auto i = std::decay_t<decltype(n)>{}; (i) != (n); ++(i))
#define REP2(i, l, r) for (auto i = (l); (i) != (r); ++(i))
#define rep(...) OVERLOAD_REP(__VA_ARGS__, REP2, REP1)(__VA_ARGS__)
#define OVERLOAD_RREP(_1, _2, _3, name, ...) name
#define RREP1(i, n) for (auto i = (n); (i) != (std::decay_t<decltype(n)>{}); --(i))
#define RREP2(i, l, r) for (auto i = (l); (i) != (r); --(i))
#define rrep(...) OVERLOAD_RREP(__VA_ARGS__, RREP2, RREP1)(__VA_ARGS__)
#define all(x) (x).begin(),(x).end()
#define rall(x) (x).rbegin(), (x).rend
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
#define in(n) insert(n)
#define mp make_pair
#define fi first
#define se second
#define print(n) cout << n << endl
#define fix(n) fixed << setprecision(n)
template <typename T>
inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false)); }
template <typename T>
inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false)); }
template<>
struct std::vector<bool>: std::basic_string<bool> {
    using std::basic_string<bool>::basic_string, std::basic_string<bool>::operator =;
    explicit vector(size_t n): vector(n, false) {}
};
const ll mod=998244353;
const ll MOD=1000000007;
const int inf = INT_MAX / 2;
const ll infl = 1LL << 60;
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
// Union-Find
struct UnionFind {
    vector<int> par, rank, siz;
    // 構造体の初期化
    UnionFind(int n) : par(n,-1), rank(n,0), siz(n,1) { }
    // 根を求める
    int root(int x) {
        if (par[x]==-1) return x; // x が根の場合は x を返す
        else return par[x] = root(par[x]); // 経路圧縮
    }
    // x と y が同じグループに属するか (= 根が一致するか)
    bool same(int x, int y) {
        return root(x)==root(y);
    }
    // x を含むグループと y を含むグループを併合する
    bool unite(int x, int y) {
        int rx = root(x), ry = root(y); // x 側と y 側の根を取得する
        if (rx==ry) return false; // すでに同じグループのときは何もしない
        // union by rank
        if (rank[rx]<rank[ry]) swap(rx, ry); // ry 側の rank が小さくなるようにする
        par[ry] = rx; // ry を rx の子とする
        if (rank[rx]==rank[ry]) rank[rx]++; // rx 側の rank を調整する
        siz[rx] += siz[ry]; // rx 側の siz を調整する
        return true;
    }
    // x を含む根付き木のサイズを求める
    int size(int x) {
        return siz[root(x)];
    }
};


#define INCLUDED_MAIN
#include __FILE__

#endif