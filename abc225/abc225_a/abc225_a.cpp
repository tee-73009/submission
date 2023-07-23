#include <bits/stdc++.h>
using namespace std;

int main() {
  string s;
  cin>>s;
  set<char> a;
  int b;
  b=s.length();
  for(int i=0;i<b;i++) a.insert(s[i]);
  if(a.size()==1)cout<<1<<endl;
  else if(a.size()==2)cout<<3<<endl;
  else cout<<6<<endl;
}