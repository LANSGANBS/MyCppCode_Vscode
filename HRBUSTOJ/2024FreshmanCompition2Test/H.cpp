#include <bits/stdc++.h>
using namespace std;
#define int long long
#define i128 __int128

i128 gcd(i128 a, i128 b)
{
    return b ? gcd(b, a % b) : a;
}

void solve()
{
    int n;
    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; ++i)
    {
        cin >> a[i];
    }
    int cs = a[0];
    for (int i = 1; i < n; ++i)
    {
        cs = gcd(cs, a[i]);
    }
    int q;
    cin >> q;
    while (q--)
    {
        int m;
        cin >> m;
        if (m % cs == 0)
        {
            cout << "YES" << endl;
        }
        else
        {
            cout << "NO" << endl;
        }
    }
}

signed main()
{
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    int tt = 1;
    cin >> tt;
    while (tt--)
    {
        solve();
    }
    return 0;
}