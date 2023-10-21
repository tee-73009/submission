#ifdef INCLUDED_MAIN

int main() {
    string s;
    cin >> s;
    int n=s.size();
    int ans = 1;
    rep(i,n){
        rep(j,i+1,n){
            int b=0;
            rep(k,j-i) if(s[i+k]!=s[j-k]) b++;
            if(b==0)ans = max(ans,j-i+1);
        }
    }
    print(ans);
}

#else
//ここからテンプレート
#include <bits/stdc++.h>
#include <atcoder/all>
#define OVERLOAD_REP(_1, _2, _3, name, ...) name
#define REP1(i, n) for (auto i = std::decay_t<decltype(n)>{}; (i) != (n); ++(i))
#define REP2(i, l, r) for (auto i = (l); (i) != (r); ++(i))
#define rep(...) OVERLOAD_REP(__VA_ARGS__, REP2, REP1)(__VA_ARGS__)
#define OVERLOAD_RREP(_1, _2, _3, name, ...) name
#define RREP1(i, n) for (auto i = (n); (i) != (std::decay_t<decltype(n)>{}); --(i))
#define RREP2(i, l, r) for (auto i = (l); (i) != (r); --(i))
#define rrep(...) OVERLOAD_RREP(__VA_ARGS__, RREP2, RREP1)(__VA_ARGS__)
#define all(x) (x).begin(),(x).end()
#define rall(x) (x).rbegin(), (x).rend()
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
#define YESNO(T) if(T){cout<<"YES"<<endl;}else{cout<<"NO"<<endl;}
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
using namespace std;
using namespace atcoder;
using ll=long long;
using ld=long double;
using ull=unsigned long long;
using vi = vector<int>;
using vl = vector<long>;
using vll = vector<long long>;
using vs = vector<string>;
using vc = vector<char>;
using vb = vector<bool>;
using vvi = vector<vi>;
using vvl = vector<vl>;
using vvc = vector<vc>;
using vvb = vector<vb>;
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
#define updiv(n,x) (n + x - 1) / x
#define print(n) cout << n << endl
#define fix(n) fixed << setprecision(n)
template <typename T>
inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false)); }
template <typename T>
inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false)); }
const ll mod=998244353;
const ll MOD=1000000007;
using mint = modint998244353;
using Mint = modint1000000007;
const long double PI=3.141592653589793;
using vm = vector<mint>;
using vvm = vector<vm>;
using vM = vector<Mint>;
using vvM = vector<vM>;
const vector<int>dx={0,1,0,-1,-1,-1,1,1};
const vector<int>dy={1,0,-1,0,-1,1,-1,1};
const int inf = INT_MAX / 2;
const ll INF= 1LL << 60;
//デバッグ
//ただの変数
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
//pair
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
//vector
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
//vector(2次元)
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
//set
#define EACH(i, s) for (__typeof__((s).begin()) i = (s).begin(); i != (s).end(); ++i)
template<class T> ostream& operator << (ostream &s, set<T> P)
{ EACH(it, P) { s << "<" << *it << "> "; } return s << endl; }
//map
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ EACH(it, P) { s << "<" << it->first << "->" << it->second << "> "; } return s << endl; }
//数学系
//繰り返し二乗法(modなし)
ll mypow(ll a, ll b) {
    if(a==0) {
        return 0;
    }
    ll res = 1;
    while(b > 0) {
        if(b & 1) res = res * a;
        a = a * a;
        b >>= 1;
    }
    return res;
}
//繰り返し二乗法(modあり)
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
// 拡張 Euclid の互除法 (a*p+b*q=1) ※p,qの宣言を忘れない　例：a=42,b=11→p=5,q=-19 (42*5+11*(-19)=1)
long long extGCD(long long a, long long b, long long& p, long long& q) {
    if (b == 0) { p = 1; q = 0; return a; }
    long long d = extGCD(b, a % b, q, p);
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
//桁和
int digsum(int n) {
    int res = 0;
    while(n > 0) {
        res += n%10;
        n /= 10;
    }
    return res;
}
//素数判定
bool isPrime(int x){
    int i;
    if(x < 2)return false;
    else if(x == 2) return true;
    if(x%2 == 0) return false;
    for(i = 3; i*i <= x; i += 2) if(x%i == 0) return false;
    return true;
}
//エラトステネスの篩（素数取得）
vector<bool> prime_table(int N) {
    vector<bool> isprime(N+1, true);
    isprime[1] = false;
    for (int p = 2; p <= N; ++p) {
        if (!isprime[p]) continue;
        for (int q = p * 2; q <= N; q += p) {
            isprime[q] = false;
        }
    }
    return isprime;
}
//素因数分解 例：360=2^3*3^2*5^1→[{2,3},{3,2},{5,1}]
vector<pll> prime_factorize(ll N) {
    vector<pll> res;
    for (ll a = 2; a * a <= N; ++a) {
        if (N % a != 0) continue;
        ll ex = 0; // 指数
 
        // 割れる限り割り続ける
        while (N % a == 0) {
            ++ex;
            N /= a;
        }
 
        // その結果を push
        res.pb({a, ex});
    }
 
    // 最後に残った数について
    if (N != 1) res.pb({N, 1});
    return res;
}
//約数全列挙
vi divisor(int n){
    vi ret;
    for(int i = 1 ; i*i <= n ; ++i){
        if(n%i == 0){
            ret.pb(i);
            if(i != 1 && i*i != n){
                ret.pb(n/i);
            }
        }
    }
    return ret;
}
//点と点の距離を返す
ld Dis(ll ax, ll ay, ll bx, ll by) {
    return sqrt(pow(ax - bx, 2) + pow(ay - by, 2));
}
//二つのベクトルから平行四辺形の面積を返す
ll he(ll x0, ll y0, ll x1, ll y1, ll x2, ll y2) {//外積を二で割る
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
}
//グラフ
//DFS
map<ll,vector<ll>>Gra;
map<ll,bool>visited;// false=not visited, true=visited
ll ans=0;
void DFS(int pos){
  ans++;
  visited[pos]=true;
  for(int i : Gra[pos]){
    if(visited[i]==false){
        DFS(i);
    }
  }
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
