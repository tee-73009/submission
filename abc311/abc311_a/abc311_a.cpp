#include <bits/stdc++.h>
using namespace std;

int main() {
  string s;
  string A="A",B="B";
  int n;
  int a=0,b=0,c=0;
  cin >> n >> s;
  for(int i=0;i<n;i++) {
    if(s[i]==A[0]) a++;
    else if(s[i]==B[0]) b++;
    else c++;
    if(a*b*c>0) {
      cout << i+1 << endl;
      break;
    }
  }
}