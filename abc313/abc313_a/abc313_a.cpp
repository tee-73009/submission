#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin>>n;
    int p[n];
    for(int i=0;i<n;i++) cin>>p[i];
    int ans=-1;
    for(int i=1;i<n;i++) ans=max(ans,p[i]-p[0]);
    cout << ans+1 << endl;
}