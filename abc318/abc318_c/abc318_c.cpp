#include <bits/stdc++.h>
#define rep(i, n) for(int i = 0; i < (int)(n); i++)
#define llrep(i, n) for(ll i = 0; i < (ll)(n); i++)
#define rrep(i, a, b) for(int i = a; i < (int)(b); i++)
#define llrrep(i, a, b) for(ll i = a; i < (ll)(b); i++)
#define all(x) (x).begin,(x).end()
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
#define YESNO(T) if(T){cout<<"YES"<<endl;}else{cout<<"NO"<<endl;}
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
using namespace std;
using ll=long long;
using ld=long double;
using ull=unsigned long long;
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
int main() {
    ll n,d,p;
    cin >> n >> d >> p;
    ll f[n];
    ll s[200009];
    rep(i,n) cin >> f[i];
    sort(f,f+n, greater<ll>());
    ll a = n/d;
    ll ans = 0;
    rep(i,n) {
        s[i+1]=s[i]+f[i]; 
    }
    rep(i,a) {
        ans+= min(s[i*d+d]-s[i*d],p);
    }
    ll t=0;
    rrep(i,a*d,n) {
        t+=f[i];
    }
    ans+=min(t,p);
    cout << ans << endl;
}