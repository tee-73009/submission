#include <bits/stdc++.h>
using namespace std;

int main() {
  string a,b;
  int x,y;
  cin>>a>>x>>y;
  b=a;
  b[x-1]=a[y-1];
  b[y-1]=a[x-1];
  cout << b << endl;
}