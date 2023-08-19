#include <bits/stdc++.h>
using namespace std;

int main(){
    int m;
    int n=1;
    int ansm=1;
    cin >> m;
    int d[m];
    for(int i=0;i<m;i++) {
        cin >> d[i];
        n+= d[i];
    }
    n/=2;
    for(int i=0;i<m;i++){
        if(d[i]>=n) {
            cout << ansm << ' '<<n<<endl;
            break;
        }
        n-=d[i];
        ansm++;
    }
}