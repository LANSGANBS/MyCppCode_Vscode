#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e5;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

int nowP;

void solve()
{
    string s;
    int x;
    cin >> s >> x;
    if (s[0] == 'P')
    {
        nowP += x;
    }
    else if (s[0] == 'B' and nowP <= x - 1)
    {
        cout << "YES" << endl;
        nowP = 0;
    }
    else
    {
        cout << "NO" << endl;
        nowP -= x;
    }
}

signed main()
{
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    int T = 1;
    cin >> T;
    while (T--)
        solve();
}