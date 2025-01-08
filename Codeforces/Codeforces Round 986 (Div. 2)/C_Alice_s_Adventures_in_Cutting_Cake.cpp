#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define i64 long long
#define all(x) x.begin(), x.end()

void solve()
{
    int n, m, v;
    cin >> n >> m >> v;
    vector<int> a(n);
    for (int i = 0; i < n; i++)
    {
        cin >> a[i];
    }
    vector<i64> s(n + 1);
    vector<int> pre(n + 1);
    vector<int> suf(n + 1);
    i64 sum = 0;
    int cnt = 0;
    for (int i = 0; i < n; i++)
    {
        s[i + 1] = s[i] + a[i];
        sum += a[i];
        if (sum >= v)
            sum = 0, cnt++;
        pre[i + 1] = cnt;
    }
    sum = 0;
    cnt = 0;
    for (int i = n - 1; ~i; i--)
    {
        sum += a[i];
        if (sum >= v)
            sum = 0, cnt++;
        suf[i] = cnt;
    }
    if (pre[n] < m)
    {
        cout << "-1\n";
        return;
    }
    i64 ans = 0;
    for (int i = n; ~i; i--)
    {
        int ned = m - suf[i];
        auto it = lower_bound(all(pre), ned) - pre.begin();
        ans = max(s[i] - s[it], ans);
    }
    cout << ans << endl;
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