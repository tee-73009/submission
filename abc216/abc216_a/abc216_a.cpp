#include <bits/stdc++.h>
using namespace std;;
int main() {
  string s;
  cin >> s;
  int n=s.length();
  int y=s[n-1]-'0';
  s=s.substr(0,n-2);
  cout << s;
  if(y<3)cout <<'-';
  if(y>6)cout <<'+';
}