#include <bits/stdc++.h>
using namespace std;
int main() {
  int n,x;
  string b="No";
  cin >> n >> x;
  vector<int>a(n);
  for(int i=0;i<n;i++){
    cin >> a[i];
  if(a[i] == x) {
    b= "Yes";
    break;
  }
  }
  cout << b << endl;
}