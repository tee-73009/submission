#include <bits/stdc++.h>
using namespace std;

int main() {
  long long x , y = 100 , a = 0;
  cin >> x;
  while (x > y) {
    a++;
    y += y / 100;
  }
  cout << a << endl;
}