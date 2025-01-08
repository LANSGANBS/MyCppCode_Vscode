#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define i64 long long

void solve()
{
    i64 n, b, c;
    cin >> n >> b >> c;
    if (n == 1)
    {
        if (c == 0)
            cout << 0 << endl;
        else
            cout << 1 << endl;
        return;
    }
    if (c >= n)
    {
        cout << n << endl;
        return;
    }
    if (b == 0)
    {
        if (c == n - 2 || c == n - 1)
        {
            cout << n - 1 << endl;
            return;
        }
        cout << "-1\n";
        return;
    }
    // b*k+c<n
    i64 k = (n - c - 1) / b + 1;
    k = min(k, n);
    cout << n - k << endl;
}

signed main()
{
    cin.tie(0);
    cout.tie(0);
    ios::sync_with_stdio(0);
    int t;
    cin >> t;
    while (t--)
        solve();
}