#include <bits/stdc++.h>
using namespace std;

int main() {
  int n ,m ,sum,y;
  sum = 0;
  cin >> n >> m;
  vector<string>c(n);
  for(int i=0;i<n;i++) {
    cin >> c[i];
  }
  vector<string>d(m);
  for(int j=0;j<m;j++) {
    cin >> d[j];
  }
  vector<int>p(m+1);
  for(int k=0;k<m+1;k++) {
    cin >> p[k];
  }
  vector<int>x(m);
  for(int l=0;l<m;l++) {
    x[l] = count(c.begin(),c.end(),d[l]);
    sum = sum + x[l] * p[l+1];
  }
  y = n - accumulate(x.begin(),x.end(),0);
  sum = sum + p[0] * y;
  cout << sum << endl;
}