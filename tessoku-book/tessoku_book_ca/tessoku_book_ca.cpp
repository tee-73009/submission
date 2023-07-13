#include <bits/stdc++.h>
using namespace std;
 
int main() {
  int a,b;
  string c;
  c = "No";
  cin >> a >> b;
  while(a <= b){
    if(100 % a == 0){
      c = "Yes";
      break;
    }
    a++;
  }
  cout << c << endl;
}