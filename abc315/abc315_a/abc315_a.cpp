#include <bits/stdc++.h>
using namespace std;

int main(){
    string s;
    char ans[109];
    int m=0;
    cin >> s;
    int n =s.length(); 
    for(int i=0;i<n;i++) {
        if(s[i]!='a'&&s[i]!='i'&&s[i]!='u'&&s[i]!='e'&&s[i]!='o') {
          ans[m]=s[i];
          m++;
        }
    }
    cout << ans << endl;
}