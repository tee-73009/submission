#include<bits/stdc++.h>
using namespace std;

int main() {
  int n,a[109];
  cin >> n;
  for(int i =0;i<n;i++) cin >>a[i];
  string ans ="No";
  for(int i=0;i<n;i++){
    for(int j=i+1;j<n;j++){
      for(int k=j+1;k<n;k++){
        if(a[i]+a[j]+a[k]==1000)ans="Yes";
      }
    }
  }
  cout << ans << endl;
}