#line 2 "library/template/template.hpp"
#include <bits/stdc++.h>
#line 3 "library/template/macro.hpp"
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
using namespace std;
#if __has_include(<atcoder/all>)
#include <atcoder/all>
using namespace atcoder;
using mint = modint998244353;
using Mint = modint1000000007;
using vm = vector<mint>;
using vvm = vector<vm>;
using vM = vector<Mint>;
using vvM = vector<vM>;
//998244353_in
istream &operator>>(istream &a, mint &b){
    long long tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//998244353_out
ostream &operator<<(ostream &a, mint &b){
    a << b.val();
    return a;
}
//1000000007_in
istream &operator>>(istream &a, Mint &b){
    long long tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//1000000007_out
ostream &operator<<(ostream &a, Mint &b){
    a << b.val();
    return a;
}
#endif
#if __has_include(<boost/multiprecision/cpp_dec_float.hpp>)
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <boost/array.hpp>
namespace mp = boost::multiprecision;
// 任意長整数型
using Bint = mp::cpp_int;
// 仮数部が32桁(10進数)の浮動小数点数型
using Real32 = mp::number<mp::cpp_dec_float<32>>;
// 仮数部が1024桁(10進数)の浮動小数点数型
using Real1024 = mp::number<mp::cpp_dec_float<1024>>;
//有理数型
using Rat = boost::rational<Bint>;
//Bint_in
istream &operator>>(istream &a, Bint &b){
    long long tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//Real32_in
istream &operator>>(istream &a, Real32 &b){
    long double tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//Real1024_in
istream &operator>>(istream &a, Real1024 &b){
    long double tmp;
    a >> tmp;
    b = tmp;
    return a;
}
#endif
#define OVERLOAD_3(_1, _2, _3, name, ...) name
#define REP1(i, n) for (auto i = std::decay_t<decltype(n)>{}; (i) != (n); ++(i))
#define REP2(i, l, r) for (auto i = (l); (i) != (r); ++(i))
#define rep(...) OVERLOAD_3(__VA_ARGS__, REP2, REP1)(__VA_ARGS__)
#define RREP1(i, n) for (auto i = (n); (i) != (std::decay_t<decltype(n)>{}); --(i))
#define RREP2(i, l, r) for (auto i = (l); (i) != (r); --(i))
#define rrep(...) OVERLOAD_3(__VA_ARGS__, RREP2, RREP1)(__VA_ARGS__)
#define rep_(i,a) rep(i,a.size())
#define all(x) std::begin(x), std::end(x)
#define rall(x) std::rbegin(x), std::rend(x)
#define yesno(T) if(T){cout<<"yes"<<endl;}else{cout<<"no"<<endl;}
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
#define YESNO(T) if(T){cout<<"YES"<<endl;}else{cout<<"NO"<<endl;}
#define exit_with(...) ({ __VA_ARGS__; return 0; })
#define break_with(...) ({ __VA_ARGS__; break; })
#define continue_with(...) ({ __VA_ARGS__; continue; })
using ll=long long;
using ld=long double;
using ull=unsigned long long;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
using vi = vector<int>;
using vll = vector<long long>;
using vs = vector<string>;
using vc = vector<char>;
using vb = vector<bool>;
using vvi = vector<vi>;
using vvll = vector<vll>;
using vvs = vector<vs>;
using vvc = vector<vc>;
using vvb = vector<vb>;
using vvvi = vector<vvi>;
using vvvll = vector<vvll>;
using vvvb = vector<vvb>;
using vpll = vector<pll>;
using mi = map<int, int>;
using mll = map<ll, ll>;
using msi = map<string, int>;
using msl = map<string, ll>;
using mmi = map<int,mi>;
using mml = map<ll,mll>;
using Graph = vector<vector<ll>>;
const string ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const string abc="abcdefghijklmnopqrstuvwxyz";
#define gcd __gcd
#define pb push_back
template<class T> using pq = priority_queue<T>;
template<class T> using pqr = priority_queue<T, vector<T>, greater<T>>;
#define lb(v, k) (lower_bound((v).begin(), (v).end(), (k)) - v.begin())
#define ub(v, k) (upper_bound((v).begin(), (v).end(), (k)) - v.begin())
#define fi first
#define se second
#define elif else if
#define updiv(n,x) (n + x - 1) / x
#define rounddiv(n,x) (ll)((double)(n)/(double)(x)+0.5)
#define fix(n) fixed << setprecision(n)
#define _print(n) cout << " " << (n) << endl
#define print_(n) cout << (n) << " "
template<class... T>
//3つ以上でも使えるmax,min
constexpr auto min(T... a){
    return min(initializer_list<common_type_t<T...>>{a...});
}
template<class... T>
constexpr auto max(T... a){
    return max(initializer_list<common_type_t<T...>>{a...});
}
template <typename T>
inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false)); }
template <typename T>
inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false)); }
const ll mod=998244353;
const ll MOD=1000000007;
const ld PI=3.141592653589793;
const vi dx4={0,-1,0,1};
const vi dy4={1,0,-1,0};
const vi dx8={1,0,-1,0,1,-1,1,-1};
const vi dy8={0,1,0,-1,1,-1,-1,1};
const vi dx6={1,0,-1,0,1,-1};
const vi dy6={0,1,0,-1,1,-1};
const int inf = INT_MAX / 2;
const ll INF= 1LL << 60;
#line 6 "library/template/inout.hpp"
#include <type_traits>
using namespace std;

template <class Char, std::enable_if_t<std::is_same_v<Char, char>, int> = 0>
inline Char in() { return getchar_unlocked(); }
template <class String, std::enable_if_t<std::is_same_v<String, std::string>, int> = 0>
inline std::string in() {
    char c; do { c = getchar_unlocked(); } while (isspace(c));
    std::string s;
    do { s.push_back(c); } while (not isspace(c = getchar_unlocked()));
    return s;
}
template <class Integer, std::enable_if_t<std::is_integral_v<Integer> and not std::is_same_v<Integer, char>, int> = 0>
inline Integer in() {
    char c; do { c = getchar_unlocked(); } while (isspace(c));
    if (std::is_signed<Integer>::value and c == '-') return -in<Integer>();
    Integer n = 0;
    do { n = n * 10 + c - '0'; } while (not isspace(c = getchar_unlocked()));
    return n;
}

template <class Char, std::enable_if_t<std::is_same_v<Char, char>, int> = 0>
inline void out(char c) { putchar_unlocked(c); }
template <class String, std::enable_if_t<std::is_same_v<String, std::string>, int> = 0>
inline void out(const std::string & s) { for (char c : s) putchar_unlocked(c); }
template <class Integer, std::enable_if_t<std::is_integral_v<Integer>, int> = 0>
inline void out(Integer n) {
    char s[20];
    int i = 0;
    if (std::is_signed<Integer>::value and n < 0) { putchar_unlocked('-'); n *= -1; }
    do { s[i ++] = n % 10; n /= 10; } while (n);
    while (i) putchar_unlocked(s[-- i] + '0');
}
#define ini(...) int __VA_ARGS__; input(__VA_ARGS__);
#define inl(...) long long __VA_ARGS__; input(__VA_ARGS__);
#define ins(...) string __VA_ARGS__; input(__VA_ARGS__);
#define inc(...) char __VA_ARGS__; input(__VA_ARGS__);
#define ing(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);name[b].pb(a);}
#define ing_on(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);}//有向
#define ing_cost(name,size,n) Graph_cost name(size); rep(i,n){inl(a,b,c);a--;b--;name[a].pb({b,c});name[b].pb({a,c});}//コスト付き
#define in1(s) for (int i = 0; i < (int)s.size();i++) input(s[i]);
#define in2(s, t) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i]);
#define in3(s, t, u) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i], u[i]);
#define in4(s, t, u, v) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i], u[i], v[i]);
void input() {}
template <typename T, class... U>
void input(T &t, U &...u) {
    t = in<T>();
    input(u...);
}
void print() { out<char>('\n'); }
template <typename T, class... U, char sep = ' '>
void print(const T &t, const U &...u) {
    out<T>(t);
    if (sizeof...(u)) out<char>(sep);
    print(u...);
}
#line 2 "library/template/func.hpp"
//XORshift
unsigned int randInt() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}
//累積和(l以上r以下→s[r+1]-s[l])
template<typename T,typename U>
void get_sum(vector<T> &a,vector<U> &sum){
    sum.resize(a.size()+1);
    rep(i,(ll)a.size()){
    sum[i+1]=sum[i]+a[i];
    }
    return;
}
//座標圧縮
// O(NlogN)
template<class T>
vector<T> compress(vector<T> &vec){
    auto vals = vec;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    for(int i = 0; i < vec.size(); i++){
        vec[i] = lower_bound(vals.begin(), vals.end(), vec[i]) - vals.begin();
    }
    return vals;
}
//繰り返し二乗法(modなし)
ull mypow(ull a, ull b) {
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
// 拡張 Euclid の互除法 (a*p+b*q=1) ※p,qの宣言を忘れない 例：a=42,b=11→p=5,q=-19 (42*5+11*(-19)=1)
ll extGCD(ll a, ll b, ll& p, ll& q) {
    if (b == 0) { p = 1; q = 0; return a; }
    ll d = extGCD(b, a % b, q, p);
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
ll div(ll a, ll b, ll Mod) {
return a * modinv(b,Mod) % Mod;
}
// 中国剰余定理
// リターン値を (r, m) とすると解は x ≡ r (mod. m)
// 解なしの場合は (0, -1) をリターン
pll ChineseRem(const vll &b, const vll &m) {
    ll r = 0, M = 1;
    for (int i = 0; i < (int)b.size(); ++i) {
        ll p, q;
        ll d = extGCD(M, m[i], p, q); // p is inv of M/d (mod. m[i]/d)
        if ((b[i] - r) % d != 0) return make_pair(0, -1);
        ll tmp = (b[i] - r) / d * p % (m[i]/d);
        r += M * tmp;
        M *= m[i]/d;
    }
    return make_pair((r%M+M)%M, M);
}
// Garner のアルゴリズム, x%MOD, LCM%MOD を求める (m は互いに素でなければならない)
// for each step, we solve "coeffs[k] * t[k] + constants[k] = b[k] (mod. m[k])"
//      coeffs[k] = m[0]m[1]...m[k-1]
//      constants[k] = t[0] + t[1]m[0] + ... + t[k-1]m[0]m[1]...m[k-2]
long long Garner(vector<long long> b, vector<long long> m, long long Mod) {
    m.push_back(Mod); // banpei
    vector<long long> coeffs((int)m.size(), 1);
    vector<long long> constants((int)m.size(), 0);
    for (int k = 0; k < (int)b.size(); ++k) {
        long long t = (((b[k] - constants[k])*modinv(coeffs[k], m[k]), m[k])+ m[k]) % m[k];
        for (int i = k+1; i < (int)m.size(); ++i) {
            (constants[i] += t * coeffs[i]) %= m[i];
            (coeffs[i] *= m[k]) %= m[i];
        }
    }
    return constants.back();
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
ll comb(ll n,ll k,ll Mod) {
return div(fact[n],(fact[k]*fact[n-k] % Mod),Mod);
}
vector<long long> fact_inv, inv; 
/*  init_nCk :二項係数のための前処理
    計算量:O(n)
*/
void init_nCk(int SIZE,int Mod) {
    fact.resize(SIZE + 5);
    fact_inv.resize(SIZE + 5);
    inv.resize(SIZE + 5);
    fact[0] = fact[1] = 1;
    fact_inv[0] = fact_inv[1] = 1;
    inv[1] = 1;
    for (int i = 2; i < SIZE + 5; i++) {
        fact[i] = fact[i - 1] * i % Mod;
        inv[i] = Mod - inv[Mod % i] * (Mod / i) % Mod;
        fact_inv[i] = fact_inv[i - 1] * inv[i] % Mod;
    }
}
/*  nCk :MODでの二項係数を求める(前処理 init_nCk が必要)
    計算量:O(1)
*/
long long nCk(int n, int k) {
    assert(!(n < k));
    assert(!(n < 0 || k < 0));
    return fact[n] * (fact_inv[k] * fact_inv[n - k] % MOD) % MOD;
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
ll isPrime(ll x) {
    ll i;
    if(x < 2)return 0;
    else if(x == 2) return 1;
    if(x%2 == 0) return 0;
    for(i = 3; i*i <= x; i += 2) if(x%i == 0) return 0;
    return 1;
}
// montgomery modint (MOD < 2^62, MOD is odd)
struct montgomerymint {
    using mint = montgomerymint;
    using u64 = uint64_t;
    using u128 = __uint128_t;
    
    // static menber
    static u64 MOD;
    static u64 INV_MOD;  // INV_MOD * MOD ≡ 1 (mod 2^64)
    static u64 T128;  // 2^128 (mod MOD)
    
    // inner value
    u64 val;
    
    // constructor
    montgomerymint() : val(0) { }
    montgomerymint(long long v) : val(reduce((u128(v) + MOD) * T128)) { }
    u64 get() const {
        u64 res = reduce(val);
        return res >= MOD ? res - MOD : res;
    }
    
    // mod getter and setter
    static u64 get_mod() { return MOD; }
    static void set_mod(u64 mod) {
        assert(mod < (1LL << 62));
        assert((mod & 1));
        MOD = mod;
        T128 = -u128(mod) % mod;
        INV_MOD = get_inv_mod();
    }
    static u64 get_inv_mod() {
        u64 res = MOD;
        for (int i = 0; i < 5; ++i) res *= 2 - MOD * res;
        return res;
    }
    static u64 reduce(const u128 &v) {
        return (v + u128(u64(v) * u64(-INV_MOD)) * MOD) >> 64;
    }
    
    // arithmetic operators
    mint operator - () const { return mint() - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        if ((val += r.val) >= 2 * MOD) val -= 2 * MOD;
        return *this;
    }
    mint& operator -= (const mint &r) {
        if ((val += 2 * MOD - r.val) >= 2 * MOD) val -= 2 * MOD;
        return *this;
    }
    mint& operator *= (const mint &r) {
        val = reduce(u128(val) * r.val);
        return *this;
    }
    mint& operator /= (const mint &r) {
        *this *= r.inv();
        return *this;
    }
    mint inv() const { return pow(MOD - 2); }
    mint pow(u128 n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }

    // other operators
    bool operator == (const mint &r) const {
        return (val >= MOD ? val - MOD : val) == (r.val >= MOD ? r.val - MOD : r.val);
    }
    bool operator != (const mint &r) const {
        return (val >= MOD ? val - MOD : val) != (r.val >= MOD ? r.val - MOD : r.val);
    }
    friend istream& operator >> (istream &is, mint &x) {
        long long t;
        is >> t;
        x = mint(t);
        return is;
    }
    friend ostream& operator << (ostream &os, const mint &x) {
        return os << x.get();
    }
    friend mint modpow(const mint &r, long long n) {
        return r.pow(n);
    }
    friend mint modinv(const mint &r) {
        return r.inv();
    }
};

typename montgomerymint::u64
montgomerymint::MOD, montgomerymint::INV_MOD, montgomerymint::T128;

// Miller-Rabin
bool MillerRabin(long long N, vector<long long> A) {
    using mint = montgomerymint;
    mint::set_mod(N);
    
    long long s = 0, d = N - 1;
    while (d % 2 == 0) {
        ++s;
        d >>= 1;
    }
    for (auto a : A) {
        if (N <= a) return true;
        mint x = mint(a).pow(d);
        if (x != 1) {
            long long t;
            for (t = 0; t < s; ++t) {
                if (x == N - 1) break;
                x *= x;
            }
            if (t == s) return false;
        }
    }
    return true;
}

bool is_prime(long long N) {
    if (N <= 1) return false;
    else if (N == 2) return true;
    else if (N % 2 == 0) return false;
    else if (N < 4759123141LL)
        return MillerRabin(N, {2, 7, 61});
    else
        return MillerRabin(N, {2, 325, 9375, 28178, 450775, 9780504, 1795265022});
}
long long pollard(long long N) {
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return N;
    
    using mint = montgomerymint;
    mint::set_mod(N);
    
    long long step = 0;
    while (true) {
        mint r = randInt();  // random r
        auto f = [&](mint x) -> mint { return x * x + r; };
        mint x = ++step, y = f(x);
        while (true) {
            long long p = gcd((y - x).get(), (uint64_t)N);
            if (p == 0 || p == N) break;
            if (p != 1) return p;
            x = f(x);
            y = f(f(y));
        }
    }
}

vector<long long> prime_factorize(long long N) {
    if (N == 1) return {};
    long long p = pollard(N);
    if (p == N) return {p};
    vector<long long> left = prime_factorize(p);
    vector<long long> right = prime_factorize(N / p);
    left.insert(left.end(), right.begin(), right.end());
    sort(left.begin(), left.end());
    return left;
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
//約数全列挙
vll divisor(ll n){
    vll ret;
    for(ll i = 1 ; i*i <= n ; ++i){
        if(n%i == 0){
            ret.pb(i);
            if(i*i != n){
                ret.pb(n/i);
            }
        }
    }
    return ret;
}
//ペル方程式
//ペル方程式の初期解を求める
pair<long long, long long> Pell(long long D, int c = 1) {
    if (D == 1) return make_pair(-1, -1);
    long long a = D, b = 1;
    while (a > b) {
        a /= 2;
        b *= 2;
    }
    a = b + 1;
    while (a > b) {
        a = b; 
        b = (b + D/b)/2;
    }
    if (a*a == D) return make_pair(-1, -1);
    
    long long a0 = a;
    bool parity = false;
    b = 1;
    long long x2 = 1, x = a, y2 = 0, y = 1, q;
    while (true) {
        b = (D - a*a) / b;
        q = (a0 + a) / b;
        a = q * b - a;
        parity = !parity;
        if (b == 1) break;
        long long tx = x, tx2 = x2, ty = y, ty2 = y2;
        x2 = tx, x = tx * q + tx2; 
        y2 = ty, y = ty * q + ty2;
    }
    long long x0 = x, y0 = y;
    if (!parity) {
        if (c == 1) return make_pair(x, y);
        else return make_pair(-1, -1);
    }
    else if (c == -1) return make_pair(x, y);
    
    long long tx = x, ty = y;
    x = x0 * tx + D * y0 * ty, y = tx * y0 + x0 * ty;
    return make_pair(x, y);
}
//初期解から一般解を求める
pair<long long, long long> PellMul(pair<long long, long long> p, pair<long long, long long> q, long long D) {
    long long f = p.first * q.first + D * p.second * q.second;
    long long s = p.first * q.second + p.second * q.first;
    return make_pair(f, s);
}
//10進数→2進数(順番が逆)
string dec2bin(ll n) {
    string r;
    while (n != 0LL) {
        r += (n % 2LL == 0LL ? "0" : "1");
        n /= 2LL;
    }
    return r;
}
//2進数→10進数
ll bin2dec(string s){
    ll r=0;
    rep(i,s.size()){
        r*=10;
        r+=(s[i]-'0');
    }
    return r;
}
//点と点の距離を返す
ld Dis(ll ax, ll ay, ll bx, ll by) {
    return sqrt(pow(ax - bx, 2) + pow(ay - by, 2));
}
//二つのベクトルから平行四辺形の面積を返す
ll he(ll x0, ll y0, ll x1, ll y1, ll x2, ll y2) {//外積を二で割る
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
}
// BIT
template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat[2];

    // [0, n)
    BIT(int n, Abel unity = 0) : UNITY_SUM(unity) {
        init(n);
    }
    void init(int n) {
        for (int iter = 0; iter < 2; ++iter)
            dat[iter].assign(n + 1, UNITY_SUM);
    }
    
    // [a, b), a and b are 0-indexed
    inline void sub_add(int p, int a, Abel x) {
        for (int i = a; i < (int)dat[p].size(); i |= i + 1)
            dat[p][i] = dat[p][i] + x;
    }
    inline void add(int a, int b, Abel x) {
        sub_add(0, a, x * (-a));
        sub_add(1, a, x);
        sub_add(0, b, x * b);
        sub_add(1, b, x * (-1));
    }
    
    // [a, b), a and b are 0-indexed
    inline Abel sub_sum(int p, int a) {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[p][i];
        return res;
    }
    inline Abel sum(int a, int b) {
        return sub_sum(0, b)
            + sub_sum(1, b) * b
            - sub_sum(0, a)
            - sub_sum(1, a) * a;
    }
    
    // debug
    void printall() {
        for (int i = 0; i < (int)dat[0].size(); ++i)
            cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};
// Segment Tree
template<class Monoid> struct SegmentTree {
    using Func = function<Monoid(Monoid, Monoid)>;

    // core member
    int N;
    Func OP;
    Monoid IDENTITY;
    
    // inner data
    int log, offset;
    vector<Monoid> dat;

    // constructor
    SegmentTree() {}
    SegmentTree(int n, const Func &op, const Monoid &identity) {
        init(n, op, identity);
    }
    SegmentTree(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init(v, op, identity);
    }
    void init(int n, const Func &op, const Monoid &identity) {
        N = n;
        OP = op;
        IDENTITY = identity;
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY);
    }
    void init(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init((int)v.size(), op, identity);
        build(v);
    }
    void pull(int k) {
        dat[k] = OP(dat[k * 2], dat[k * 2 + 1]);
    }
    void build(const vector<Monoid> &v) {
        assert(N == (int)v.size());
        for (int i = 0; i < N; ++i) dat[i + offset] = v[i];
        for (int k = offset - 1; k > 0; --k) pull(k);
    }
    int size() const {
        return N;
    }
    Monoid operator [] (int i) const {
        return dat[i + offset];
    }
    
    // update A[i], i is 0-indexed, O(log N)
    void set(int i, const Monoid &v) {
        assert(0 <= i && i < N);
        int k = i + offset;
        dat[k] = v;
        while (k >>= 1) pull(k);
    }
    
    // get [l, r), l and r are 0-indexed, O(log N)
    Monoid prod(int l, int r) {
        assert(0 <= l && l <= r && r <= N);
        Monoid val_left = IDENTITY, val_right = IDENTITY;
        l += offset, r += offset;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) val_left = OP(val_left, dat[l++]);
            if (r & 1) val_right = OP(dat[--r], val_right);
        }
        return OP(val_left, val_right);
    }
    Monoid all_prod() {
        return dat[1];
    }
    
    // get max r such that f(v) = True (v = prod(l, r)), O(log N)
    // f(IDENTITY) need to be True
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == N) return N;
        l += offset;
        Monoid sum = IDENTITY;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(OP(sum, dat[l]))) {
                while (l < offset) {
                    l = l * 2;
                    if (f(OP(sum, dat[l]))) {
                        sum = OP(sum, dat[l]);
                        ++l;
                    }
                }
                return l - offset;
            }
            sum = OP(sum, dat[l]);
            ++l;
        } while ((l & -l) != l);  // stop if l = 2^e
        return N;
    }

    // get min l that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int min_left(const function<bool(Monoid)> f, int r = -1) {
        if (r == 0) return 0;
        if (r == -1) r = N;
        r += offset;
        Monoid sum = IDENTITY;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(OP(dat[r], sum))) {
                while (r < offset) {
                    r = r * 2 + 1;
                    if (f(OP(dat[r], sum))) {
                        sum = OP(dat[r], sum);
                        --r;
                    }
                }
                return r + 1 - offset;
            }
            sum = OP(dat[r], sum);
        } while ((r & -r) != r);
        return 0;
    }
    
    // debug
    friend ostream& operator << (ostream &s, const SegmentTree &seg) {
        for (int i = 0; i < (int)seg.size(); ++i) {
            s << seg[i];
            if (i != (int)seg.size() - 1) s << " ";
        }
        return s;
    }
};
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
//DFS
map<ll,vector<ll>>Gra;
map<ll,bool>visited;// false=not visited, true=visited
ll ans=0;
void DFS(ll pos){
    ans++;
    visited[pos]=true;
    for(ll i : Gra[pos]){
        if(visited[i]==false){
            DFS(i);
        }
    }
}
//BFS
Graph BFS(ll H, ll W, const vector<string> &G, pll s) {
    vector<vector<ll>> dist(H, vector<ll>(W, -1));  //すべての頂点を未訪問に初期化
    queue<pll> que;
    //初期条件 (頂点sを初期頂点とする)
    dist[s.first][s.second] = 0;
    que.push(s);  // sを探索済み頂点に
    // BFS開始
    while (!que.empty()) {
        pll v = que.front();
        que.pop();
        //頂点vからたどれる頂点を全て調べる
        for (ll i = 0; i < 4; i++) {
            ll X = dx4[i] + v.first;
            ll Y = dy4[i] + v.second;
            if (X < 0 || X >= H || Y < 0 || Y >= W) continue;
            //すでに発見済みの頂点は探索しない
            if (dist[X][Y] != -1 || G[X][Y] == '#') continue;
            //新たな未探索点xについて距離情報を更新してキューに挿入
            dist[X][Y] = dist[v.first][v.second] + 1;
            que.push(make_pair(X, Y));
        }
    }
    return dist;
}
struct edge {
    long long to;
    long long cost;
};
using Graph_cost = vector<vector<edge>>;
void dijkstra(const Graph_cost &G, int s, vll &dis) {
    ll n = G.size();
    pqr<pll> pq;  // 「仮の最短距離, 頂点」が小さい順に並ぶ
    dis.resize(n, INF);
    dis[s] = 0;
    pq.emplace(0, s);
    while (!pq.empty()) {
        pll p = pq.top();
        pq.pop();
        ll v = p.second;
        if (dis[v] != -1) {  // 最短距離で無ければ無視
            for (auto &e : G[v]) {
                long long dp = dis[v] + e.cost;
                if (dis[e.to] > dp) {  // 最短距離候補なら priority_queue に追加
                    dis[e.to] = dp;
                    pq.emplace(dp, e.to);
                }
            }
        }
    }
}
// ローリングハッシュ
struct RollingHash {
    static const long long base1 = 1007, base2 = 2009;
    static const long long mod1 = 1000000007, mod2 = 1000000009;
    vector<long long> hash1, hash2, power1, power2;

    // construct
    RollingHash(const string &S) {
        long long n = (long long)S.size();
        hash1.assign(n+1, 0), hash2.assign(n+1, 0);
        power1.assign(n+1, 1), power2.assign(n+1, 1);
        for (long long i = 0; i < n; ++i) {
            hash1[i+1] = (hash1[i] * base1 + S[i]) % mod1;
            hash2[i+1] = (hash2[i] * base2 + S[i]) % mod2;
            power1[i+1] = (power1[i] * base1) % mod1;
            power2[i+1] = (power2[i] * base2) % mod2;
        }
    }
    
    // get hash value of S[left:right]
    inline long long get(long long l, long long r) const {
        long long res1 = hash1[r] - hash1[l] * power1[r-l] % mod1;
        if (res1 < 0) res1 += mod1;
        long long res2 = hash2[r] - hash2[l] * power2[r-l] % mod2;
        if (res2 < 0) res2 += mod2;
        return res1 * mod2 + res2;
    }

    // get hash value of S
    inline long long get() const {
        return hash1.back() * mod2 + hash2.back();
    }

    // get lcp of S[a:] and S[b:]
    inline long long getLCP(long long a, long long b) const {
        long long len = min((long long)hash1.size()-a, (long long)hash1.size()-b);
        long long low = 0, high = len;
        while (high - low > 1) {
            long long mid = (low + high) >> 1;
            if (get(a, a+mid) != get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }

    // get lcp of S[a:] and T[b:]
    inline long long getLCP(const RollingHash &T, long long a, long long b) const {
        long long len = min((long long)hash1.size()-a, (long long)hash1.size()-b);
        long long low = 0, high = len;
        while (high - low > 1) {
            long long mid = (low + high) >> 1;
            if (get(a, a+mid) != T.get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }
};
#line 6 "library/template/template.hpp"
using namespace std;
#line 2 "code.cpp"
#if __has_include(<atcoder/all>)
#include<atcoder/all>
using namespace atcoder;
#endif
signed main() {
    inl(h,w,n);
    vll a(n);
    in1(a);
    sort(all(a));
    ll x=a[0];
    if(x!=0){
    rep(i,n) a[i]-=x;
    h=h>>(x);
    w=w>>(x);
    }
    sort(rall(a));
    ll TMP=max(h,w);
    ll tmp=min(h,w);
    h=TMP;
    w=tmp;
    string H=dec2bin(h);
    string W=dec2bin(w);
    ll hs=H.size();
    ll ws=W.size();
    vll A(a[0]+1);
    rep(i,n){
        A[a[i]]++;
    }
    bool can=true;
    for(ll i=a[0];i>=0;i--){
        ll x=h>>(i);
        ll y=w>>(i);
        ll z=x*y;
        if(z<A[i]) can=false;
        if(i!=0) A[i-1]+=(4*A[i]);
    }
    YesNo(can);
}
