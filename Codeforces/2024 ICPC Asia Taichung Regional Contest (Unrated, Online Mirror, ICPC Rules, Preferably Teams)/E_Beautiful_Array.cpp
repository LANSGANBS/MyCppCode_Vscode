#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e5;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

void solve()
{
    int a, b;
    cin >> a >> b;

    if (a == b)
    {
        cout << 1 << endl;
        cout << a << endl;
    }
    else
    {
        cout << 3 << endl;
        if (a < b)
        {
            int x = 3 * a - 2 * b;
            cout << x << ' ' << b << ' ' << b << endl;
        }
        else
        {
            int y = 3 * a - 2 * b;
            cout << b << ' ' << b << ' ' << y << endl;
        }
    }
}

signed main()
{
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    int T = 1;
    // cin >> T;
    while (T--)
        solve();
}