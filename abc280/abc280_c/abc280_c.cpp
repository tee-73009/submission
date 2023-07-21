#include <bits/stdc++.h>
using namespace std;

int main() {
  string s,t;
  cin >> s >> t;
  for(int i=0;i<500009;i++) {
    if(s[i]!=t[i]) {
      cout << i+1 << endl;
      break;
    }
  }
}