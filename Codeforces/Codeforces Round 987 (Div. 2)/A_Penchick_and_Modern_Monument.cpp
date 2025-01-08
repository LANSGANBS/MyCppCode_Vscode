#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

#define all(x) x.begin(), x.end()

void solve()
{
    int n;
    cin >> n;
    vector<int> dp(n, 1);
    vector<int> a(n);
    for (int i = 0; i < n; i++)
        cin >> a[i];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (a[j] <= a[i])
                dp[i] = max(dp[i], dp[j] + 1);
        }
    }
    cout << n - *max_element(all(dp)) << endl;
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