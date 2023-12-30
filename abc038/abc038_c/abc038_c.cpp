#ifdef INCLUDED_MAIN

signed main() {
    inl(n);
    vll a(n);
    IN(a[0]);
    ll ans=0;
    ll x=1;
    rep(i,n-1){
        IN(a[i+1]);
        if(a[i]<a[i+1]) x++;
        else{
            ans+=(x*(x+1)/2);
            x=1;
        }
    }
    ans+=(x*(x+1)/2);
    print(ans);
}


#else
#pragma region template //テンプレート
#pragma region macro //マクロ
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
#define rep_(i,a) rep(i,a.size())
#define all(x) (x).begin(),(x).end()
#define rall(x) (x).rbegin(), (x).rend()
#define yesno(T) if(T){cout<<"yes"<<endl;}else{cout<<"no"<<endl;}
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
#define YESNO(T) if(T){cout<<"YES"<<endl;}else{cout<<"NO"<<endl;}
#pragma endregion
#pragma region boost::multiprecision //多倍長整数
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
#pragma endregion
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
using vvll = vector<vll>;
using vvs = vector<vs>;
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
#define pq(T) priority_queue<T>
#define pqr(T) priority_queue<T, vector<T>, greater<T>>
#define lb(v, k) (lower_bound((v).begin(), (v).end(), (k)) - v.begin())
#define ub(v, k) (upper_bound((v).begin(), (v).end(), (k)) - v.begin())
#define fi first
#define se second
#define elif else if
#define updiv(n,x) (n + x - 1) / x
#define rounddiv(n,x) (ll)((double)(n)/(double)(x)+0.5)
#define print(n) cout << (n) << endl
#define _print(n) cout << " " << (n) << endl
#define print_(n) cout << (n) << " "
#define fix(n) fixed << setprecision(n)
template <typename T>
inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false)); }
template <typename T>
inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false)); }
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
// index が条件を満たすかどうか
bool isOK(vll a,ll index,ll key) {
    if (a[index] >= key) return true;
    else return false;
}
// 汎用的な二分探索のテンプレ
ll binary_search(vll a,ll key) {
    ll ng = -1; //「index = 0」が条件を満たすこともあるので、初期値は -1
    ll ok = (ll)a.size(); // 「index = a.size()-1」が条件を満たさないこともあるので、初期値は a.size()

    /* ok と ng のどちらが大きいかわからないことを考慮 */
    while (abs(ok - ng) > 1) {
        ll mid = (ok + ng) / 2;

        if (isOK(a,mid, key)) ok = mid;
        else ng = mid;
    }
    return ok;
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


const ll mod=998244353;
const ll MOD=1000000007;
using mint = modint998244353;
using Mint = modint1000000007;
const ld PI=3.141592653589793;
using vm = vector<mint>;
using vvm = vector<vm>;
using vM = vector<Mint>;
using vvM = vector<vM>;
const vi dx4={1,0,-1,0};
const vi dy4={0,1,0,-1};
const vi dx8={1,0,-1,0,1,-1,1,-1};
const vi dy8={0,1,0,-1,1,-1,-1,1};
const vi dx6={1,0,-1,0,1,-1};
const vi dy6={0,1,0,-1,1,-1};
const int inf = INT_MAX / 2;
const ll INF= 1LL << 60;
#pragma endregion
#pragma region inout //入力・出力
#define ini(...) int __VA_ARGS__; IN(__VA_ARGS__)
#define inl(...) long long __VA_ARGS__; IN(__VA_ARGS__)
#define ins(...) string __VA_ARGS__; IN(__VA_ARGS__)
#define inc(...) char __VA_ARGS__; IN(__VA_ARGS__)
#define ing(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);name[b].pb(a);}
#define ing_on(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);}//有向
#define ing_cost(name,size,n) Graph_cost name(size); rep(i,n){inl(a,b,c);a--;b--;name[a].pb({b,c});name[b].pb({a,c});}//コスト付き
#define in1(s) for (int i = 0; i < (int)s.size();i++) IN(s[i]);
#define in2(s, t) for (int i = 0; i < (int)s.size(); i++) IN(s[i], t[i]);
#define in3(s, t, u) for (int i = 0; i < (int)s.size(); i++) IN(s[i], t[i], u[i]);
#define in4(s, t, u, v) for (int i = 0; i < (int)s.size(); i++) IN(s[i], t[i], u[i], v[i]);
//pair_out
template <typename T, typename U>
ostream &operator<<(ostream &os, const pair<T, U> &p) {
  os << p.first << " " << p.second;
  return os;
}
//pair_in
template <typename T, typename U>
istream &operator>>(istream &is, pair<T, U> &p) {
  is >> p.first >> p.second;
  return is;
}
//vector_out
template <typename T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  int s = (int)v.size();
  for (int i = 0; i < s; i++) os << (i ? " " : "") << v[i];
  return os;
}
//vector_in
template <typename T>
istream &operator>>(istream &is, vector<T> &v) {
  for (auto &x : v) is >> x;
  return is;
}
//でかい数の高速化
istream &operator>>(istream &is, __int128_t &x) {
  string S;
  is >> S;
  x = 0;
  int flag = 0;
  for (auto &c : S) {
    if (c == '-') {
      flag = true;
      continue;
    }
    x *= 10;
    x += c - '0';
  }
  if (flag) x = -x;
  return is;
}
istream &operator>>(istream &is, __uint128_t &x) {
  string S;
  is >> S;
  x = 0;
  for (auto &c : S) {
    x *= 10;
    x += c - '0';
  }
  return is;
}
ostream &operator<<(ostream &os, __int128_t x) {
  if (x == 0) return os << 0;
  if (x < 0) os << '-', x = -x;
  string S;
  while (x) S.push_back('0' + x % 10), x /= 10;
  reverse(begin(S), end(S));
  return os << S;
}
ostream &operator<<(ostream &os, __uint128_t x) {
  if (x == 0) return os << 0;
  string S;
  while (x) S.push_back('0' + x % 10), x /= 10;
  reverse(begin(S), end(S));
  return os << S;
}
//vector<vector>_out
template <typename T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v)
{
    for (int i = 0; i < (int)v.size(); i++)
    {
        os << v[i] << endl;
    }
    return os;
}
//vector<vector<vector>>_out
template <typename T>
ostream &operator<<(ostream &os, const vector<vector<vector<T>>> &v)
{
    for (int i = 0; i < (int)v.size(); i++)
    {
        os << "i = " << i << endl;
        os << v[i];
    }
    return os;
}
//map_out
/*
template <typename T, typename S>
ostream &operator<<(ostream &os, const map<T, S> &mp)
{
    for ((auto &[key, val]) : mp)
    {
        os << key << ":" << val << " ";
    }
    return os;
}
*/
//set_out
template <typename T>
ostream &operator<<(ostream &os, const set<T> &st)
{
    auto itr = st.begin();
    for (int i = 0; i < (int)st.size(); i++)
    {
        os << *itr << (i + 1 != (int)st.size() ? " " : "");
        itr++;
    }
    return os;
}
//multiset_out
template <typename T>
ostream &operator<<(ostream &os, const multiset<T> &st)
{
    auto itr = st.begin();
    for (int i = 0; i < (int)st.size(); i++)
    {
        os << *itr << (i + 1 != (int)st.size() ? " " : "");
        itr++;
    }
    return os;
}
//queue_out
template <typename T>
ostream &operator<<(ostream &os, queue<T> q)
{
    while (q.size())
    {
        os << q.front() << " ";
        q.pop();
    }
    return os;
}
//deque_out
template <typename T>
ostream &operator<<(ostream &os, deque<T> q)
{
    while (q.size())
    {
        os << q.front() << " ";
        q.pop_front();
    }
    return os;
}
//stack_out
template <typename T>
ostream &operator<<(ostream &os, stack<T> st)
{
    while (st.size())
    {
        os << st.top() << " ";
        st.pop();
    }
    return os;
}
//priority_queue_out
template <class T, class Container, class Compare>
ostream &operator<<(ostream &os, priority_queue<T, Container, Compare> pq)
{
    while (pq.size())
    {
        os << pq.top() << " ";
        pq.pop();
    }
    return os;
}
//998244353_in
istream &operator>>(istream &a, mint &b){
    ll tmp;
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
    ll tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//1000000007_out
ostream &operator<<(ostream &a, Mint &b){
    a << b.val();
    return a;
}
void IN() {}
template <typename T, class... U>
void IN(T &t, U &...u) {
  cin >> t;
  IN(u...);
}

void out() { cout << "\n"; }
template <typename T, class... U, char sep = ' '>
void out(const T &t, const U &...u) {
  cout << t;
  if (sizeof...(u)) cout << sep;
  out(u...);
}
#pragma endregion
#pragma region debug //デバッグ
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
#pragma endregion
#pragma region math //数学系
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
// 拡張 Euclid の互除法 (a*p+b*q=1) ※p,qの宣言を忘れない　例：a=42,b=11→p=5,q=-19 (42*5+11*(-19)=1)
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
long long Garner(vector<long long> b, vector<long long> m, long long MOD) {
    m.push_back(MOD); // banpei
    vector<long long> coeffs((int)m.size(), 1);
    vector<long long> constants((int)m.size(), 0);
    for (int k = 0; k < (int)b.size(); ++k) {
        long long t = (((b[k] - constants[k]) * modinv(coeffs[k], m[k]), m[k])+ m[k]) % m[k];
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
ll comb(ll n,ll k,ll mod) {
return div(fact[n],(fact[k]*fact[n-k] % mod),mod);
}
vector<long long> fact_inv, inv; 
/*  init_nCk :二項係数のための前処理
    計算量:O(n)
*/
void init_nCk(int SIZE,int MOD) {
    fact.resize(SIZE + 5);
    fact_inv.resize(SIZE + 5);
    inv.resize(SIZE + 5);
    fact[0] = fact[1] = 1;
    fact_inv[0] = fact_inv[1] = 1;
    inv[1] = 1;
    for (int i = 2; i < SIZE + 5; i++) {
        fact[i] = fact[i - 1] * i % MOD;
        inv[i] = MOD - inv[MOD % i] * (MOD / i) % MOD;
        fact_inv[i] = fact_inv[i - 1] * inv[i] % MOD;
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
ull isPrime(ull x) {
    ll i;
    if(x < 2)return 0;
    else if(x == 2) return 1;
    if(x%2 == 0) return 0;
    for(i = 3; i*i <= x; i += 2) if(x%i == 0) return i;
    return 1;
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
//10進数→2進数
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
#pragma endregion
#pragma region data structure
template<typename T> class BIT {
private:
	int n;
	vector<T> bit;
public:
	// 0_indexed で i 番目の要素に x を加える
	void add(int i, T x){
		i++;
		while(i < n){
			bit[i] += x, i += i & -i;
		}
	}
	// 0_indexed で [0,i] の要素の和(両閉区間)
	T sum_sub(int i){
		i++;
		T s = 0;
		while(i > 0){
			s += bit[i], i -= i & -i;
		}
		return s;
	}
    T sum(int i,int j){
        return sum_sub(j)-sum_sub(i-1);
    }
    //[0,i] の要素の和>=xとなる最小のiを求める(任意のkでbit[k]>=0が必要)
    int lower_bound(T x){
        if(x<=0) return 0;
        else {
            int i=0,r=1;
            while(r<n) r=r<<1;
            for(int len=r;len>0;len=len>>1) {
                if(i+len<n && bit[i+len]<x){
                    x-=bit[i+len];
                    i+=len;
                }
            }
        }
    }
	BIT(){}
	//初期値がすべて0の場合
	BIT(int sz) : n(sz+1), bit(n, 0){}
	BIT(const vector<T>& v) : n((int)v.size()+1), bit(n, 0){
		for(int i = 0; i < n-1; i++){
			add(i,v[i]);
		}
	}
	void printall(){
		for(int i = 0; i < n-1; i++){
			cout << sum(i) - sum(i-1) << " ";
		}
		cout << "\n";
	}
	//-1スタート
	void print_sum(){
		for(int i = 0; i < n; i++){
			cout << sum(i-1) << " ";
		}
		cout << "\n";
	}
};
 
// u を昇順にソートするのに必要な交換回数(転倒数) (u は {0,..., n-1} からなる重複を許した長さ n の数列)
ll inv_count(const vector<long long>& u,const int Max)
{
    int n=u.size();
    int m=Max;
	BIT<long long> bt(m);
	long long ans = 0;
	for(int i = 0; i < n; i++){
		ans += i - bt.sum_sub(u[i]);
		bt.add(u[i], 1);
	}
	return ans;
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
#pragma endregion
#pragma region graph//グラフ
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
    pqr(pll) pq;  // 「仮の最短距離, 頂点」が小さい順に並ぶ
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
#pragma endregion
#pragma region string//文字列
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


#pragma endregion

#define INCLUDED_MAIN
#include __FILE__

#endif