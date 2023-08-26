#include <bits/stdc++.h>
using namespace std;

int main()
{
    int n,h,x;
    cin>>n>>h>>x;
    int a=x-h;
    int p[n];
    for(int i=0;i<n;i++) {
        cin >> p[i];
        if(p[i]>=a) {
          cout << i+1 << endl;
          break;
        }
    }
}