#include <bits/stdc++.h>
using namespace std;

int main()
{
    int n;
    cin >> n;
    int a[n];
    for(int i=0;i<n;i++) {
        cin >> a[i];
    }
    priority_queue<int ,vector<int>,greater<int>> b;
    for(int i=0;i<n;i++) b.push(a[i]);
    int p,q;
    p=b.top();
    b.pop();
    for(int i=0;i<n;i++) {
        q=b.top();
        if(q-p!=1){
            cout << p+1 << endl;
            break;
        }
        b.pop();
        p=q;
    }
}