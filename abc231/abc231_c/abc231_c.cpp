#include <bits/stdc++.h>
using namespace std;

int main()
{
    int n,q;
    cin>>n>>q;
    int a[n];
    for(int i=0;i<n;i++)cin>>a[i];
    int x[q];
    for(int i=0;i<q;i++)cin>>x[i];
    sort(a,a + n);
    for(int i=0;i<q;i++)cout << n - (lower_bound(a,a+n,x[i])-a)<<endl;
}