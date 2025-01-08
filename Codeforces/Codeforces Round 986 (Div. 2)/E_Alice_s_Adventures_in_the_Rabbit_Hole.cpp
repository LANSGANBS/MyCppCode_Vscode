#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long
#define all(x) begin(x), end(x)
const int N = 1e6 + 10;
const int mod = 998244353;
#define i64 long long

i64 qp(i64 a, i64 b)
{
    i64 res = 1;
    a %= mod;
    while (b)
    {
        if (b & 1)
            res = res * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return res;
}

i64 dp[N];
void init()
{
    // dp[i]=dp[i-1]*0.5+dp[i+1]*0.5
    i64 inv2 = qp(2, mod - 2);
    dp[1] = inv2;
    for (int j = 2; j < N; j++)
    {
        i64 v = (1 - dp[j - 1] * inv2 % mod + mod) % mod;
        dp[j] = inv2 * qp(v, mod - 2) % mod;
    }
}
void solve()
{
    int n;
    cin >> n;
    vector<vector<int>> e(n);
    for (int i = 1; i < n; i++)
    {
        int u, v;
        cin >> u >> v;
        u--, v--;
        e[u].push_back(v);
        e[v].push_back(u);
    }
    vector<int> mn(n);
    function<void(int x, int fa)> dfs = [&](int x, int fa)
    {
        mn[x] = n;
        bool flag = false;
        for (auto y : e[x])
        {
            if (y == fa)
                continue;
            dfs(y, x);
            flag = true;
            mn[x] = min(mn[y] + 1, mn[x]);
        }
        if (!flag)
            mn[x] = 0;
    };
    dfs(0, -1);
    vector<i64> ans(n);
    function<void(int x, int fa)> dfs2 = [&](int x, int fa)
    {
        for (auto y : e[x])
        {
            if (y == fa)
                continue;
            // dp[j]=dp[j+1]*
            ans[y] = dp[mn[y]] * ans[x] % mod;
            dfs2(y, x);
        }
    };
    ans[0] = 1;
    dfs2(0, -1);
    for (int i = 0; i < n; i++)
    {
        cout << ans[i] << " \n"[i == n - 1];
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