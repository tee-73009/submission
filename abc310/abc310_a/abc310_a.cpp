#include <bits/stdc++.h>
using namespace std;
int main() {
  int n,p,q;
  cin>>n>>p>>q;
  int ans=p;
  int d[n];
  for(int i=0;i<n;i++) cin>>d[i];
  for(int i=0;i<n;i++) ans=min(ans,q+d[i]);
  cout << ans << endl;
}