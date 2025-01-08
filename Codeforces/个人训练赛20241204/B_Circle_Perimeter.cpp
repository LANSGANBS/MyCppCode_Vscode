#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

void solve()
{
    int x, ans = 0;
    cin >> x;
    int r = x, l = x;
    for (int i = 0; i <= x; i++)
    {
        while (i * i + r * r >= (x + 1) * (x + 1))
        {
            r--;
        }
        while (l > 0 and i * i + (l - 1) * (l - 1) >= (x) * (x))
        {
            l--;
        }
        if (i * i + r * r >= (x) * (x) and i * i + l * l < (x + 1) * (x + 1))
        {
            ans += (r - l + 1);
        }
    }
    cout << ans * 4 - 4 << endl;
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