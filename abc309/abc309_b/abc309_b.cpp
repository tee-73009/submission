#include <bits/stdc++.h>
using namespace std;
int main() {
  int n;
  cin >> n;
  vector<string> x(n);
  vector<vector<int>> a(n, vector<int>(n));
  vector<vector<int>> b(n, vector<int>(n));
  for (int i=0;i<n;i++) {
    cin >> x.at(i) ;
  }
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      a.at(i).at(j) = (x.at(i))[j];
      a.at(i).at(j) = a.at(i).at(j) - 48;
      b.at(i).at(j) = a.at(i).at(j);
    }
  }
  for (int i=0;i<n-1;i++) {
    b.at(0).at(i+1)=a.at(0).at(i);
    b.at(i+1).at(n-1)=a.at(i).at(n-1);
    b.at(n-1).at(i)=a.at(n-1).at(i+1);
    b.at(i).at(0)=a.at(i+1).at(0);
  }
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      cout << b.at(i).at(j);
    }
    cout << endl;
  }
}
    