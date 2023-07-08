#include <bits/stdc++.h>
using namespace std;
int main() {
  int a,b;
  cin >> a >> b;
  if (b-a != 1) {
    cout << "No" << endl;
  } else if (a == 3 || a == 6) {
    cout << "No" << endl;
  } else {
    cout << "Yes" << endl;
  }
}