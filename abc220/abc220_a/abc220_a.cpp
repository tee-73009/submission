#include <bits/stdc++.h>
using namespace std;

int main() {
  int a,b,c,d=0;
  cin>>a>>b>>c;
  for(int i=a;i<=b;i++) {
    if(i%c==0) {
      cout <<i<<endl;
      d++;
      break;
    }
  }
  if(d==0) cout << -1<< endl;
}