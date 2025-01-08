#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

void solve()
{
    int n;
    cin >> n;
    if (n % 2 == 0)
    {
        vector<int> ans(n);
        int cnt = 0;
        for (int i = 0; i < n; i += 2)
        {
            ans[i] = ans[i + 1] = ++cnt;
        }
        for (int i = 0; i < n; i++)
        {
            cout << ans[i] << " \n"[i == n - 1];
        }
        return;
    }
    else if (n >= 27)
    {
        vector<int> ans(n);
        int cnt = 0;
        ans[0] = ans[9] = ans[25] = ++cnt;
        ans[22] = ans[26] = ++cnt;
        ans[23] = ans[24] = ++cnt;
        for (int i = 1; i < 9; i += 2)
        {
            ans[i] = ans[i + 1] = ++cnt;
        }
        for (int i = 10; i < 22; i += 2)
        {
            ans[i] = ans[i + 1] = ++cnt;
        }
        for (int i = 27; i < n; i += 2)
        {
            ans[i] = ans[i + 1] = ++cnt;
        }
        for (int i = 0; i < n; i++)
        {
            cout << ans[i] << " \n"[i == n - 1];
        }
    }
    else
    {
        cout << "-1\n";
        return;
    }
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