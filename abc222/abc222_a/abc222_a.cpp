#include<bits/stdc++.h>
using namespace std;

int main() {
  int n;
  int x=0;
  cin >> n;
  if (n==0)cout<<"0000"<<endl;
  else{
  	for(int i=n;i<10000;i*=10)x++;
  	for(int i=0;i<x-1;i++)cout<<0;
  	cout<<n<<endl;
  }
}