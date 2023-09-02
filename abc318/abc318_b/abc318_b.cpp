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
    set<int> p;
    int n;
    cin>>n;
    int a[n],b[n],c[n],d[n];
    rep(i,n) cin >> a[i] >> b[i] >> c[i] >> d[i];
    rep(i,n){
        rrep(j,a[i],b[i]){
            rrep(k,c[i],d[i]){
                p.insert(1000*j+k);
            }
        }
    }
    cout << p.size() << endl;
}