#include<bits/stdc++.h>
using namespace std;

int main() {
  string a,ans;
  cin >> a;
  string b ="Hello,World!";
  ans="AC";
  if(a.length()!=12) {
    ans="WA";
  }
  for(int i = 0;i<12;i++) {
    if(a[i]!=b[i]) {
      ans="WA";
    }
  }
  cout << ans << endl;
}