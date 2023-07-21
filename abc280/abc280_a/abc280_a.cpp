#include <bits/stdc++.h>
using namespace std;

int main() {
  int h,w;
  string s;
  cin >> h >> w;
  int ans = 0;
  for(int i=0;i<h;i++) {
    cin >> s;
    for(int j=0;j<w;j++) {
      if(s[j]=='#') ans++;
    }
  }
  cout << ans << endl;
}