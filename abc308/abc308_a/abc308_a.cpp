#include <bits/stdc++.h>
using namespace std;

int main() {
  int s1,s2,s3,s4,s5,s6,s7,s8;
  cin >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8;
  if (s1 % 25 != 0 || s2 % 25 != 0 || s3 % 25 != 0 || s4 % 25 != 0 || s5 % 25 != 0 || s6 % 25 != 0 || s7 % 25 != 0 || s8 % 25 != 0) {
    cout << "No" << endl;
  }else if (100 > s1 || s1 > s2 || s2 > s3 || s3 > s4 || s4 > s5 || s5 > s6 || s6 > s7 || s7 > s8 || s8 > 675) {
    cout << "No" << endl;
  }else {
    cout << "Yes" << endl;
  }
}