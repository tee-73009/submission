#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#include <immintrin.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <cfenv>
#include <cfloat>
#include <chrono>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <streambuf>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
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
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
using namespace std;
signed main() {
    ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    short n,m;
    cin >> n >> m;
    short a[n][m];
    rep(i,n) rep(j,m) cin >> a[i][j];
    int ans=0;
    for(short i=0;i<n;i++){
        for(short j=i+1;j<n;j++){
            int b=0;
            for(short k=0;k<m;k++){
                b^=(a[i][k]==a[j][k]);
            }
            ans+=b;
        }
    }
    cout << ans << endl;
}