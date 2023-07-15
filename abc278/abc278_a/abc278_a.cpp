#include <bits/stdc++.h>
using namespace std;

int main() {
  int n,k,a[109];
  cin >> n >> k;
  for(int i=0;i<n;i++)cin>>a[i];
  if(n>k) {
    for(int i=k;i<n;i++)cout<<a[i]<<" ";
    for(int i=0;i<k-1;i++)cout<<0<<" ";
  	cout << 0 << endl;
}else{
    for(int i=0;i<n-1;i++)cout<<0<<" ";
  	cout << 0 << endl;
  }
}