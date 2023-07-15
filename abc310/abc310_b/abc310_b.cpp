#include <bits/stdc++.h>
using namespace std;

int main() {
  int n,m,x,y;
  string ans ="No";
  cin >> n >> m;
  vector<int> p(n);
  vector<int> c(n);
  vector<vector<int>> f(109,vector<int>(109));
  for(int i=0;i<n;i++) {
    cin>>p[i]>>c[i];
    for(int j=0;j<c[i];j++) cin>>f[i][j];
  }
  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++) {
      set<int> a;
      for(int k=0;k<c[j];k++) {
        for(int l=0;l<c[i];l++) {
          a.insert(f[i][l]);
          x=a.size();
        }
        a.insert(f[j][k]);
        y=a.size();
        if(x!=y) goto next;
      }
      if(x==y&&c[i]>c[j]) {
        if (p[i]<=p[j]) ans = "Yes";
      } else if(x==y&&c[i]==c[j]) {
        if (p[i]!=p[j]) ans = "Yes";
      }
      next:;
    }
  }
  cout << ans << endl;
}
      