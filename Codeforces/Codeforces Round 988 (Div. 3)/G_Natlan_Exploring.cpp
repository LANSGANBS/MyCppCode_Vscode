#include <bits/stdc++.h>
using namespace std;
using i64 = long long;
const int N = 1e6 + 10;
const int mod = 998244353;
#define endl '\n'
bool flag = false;
int cnt = 0;
#define all(x) x.begin(), x.end()
// int xx[] = { 1,0,-1,0 };
// int yy[] = { 0,1,0,-1 };
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

void solve()
{
    int n;
    cin >> n;
    vector<int> a(n + 1);
    for (int i = 1; i <= n; ++i)
    {
        cin >> a[i];
    }

    vector<int> dp(n + 1, 0);
    dp[1] = 1;

    for (int i = 1; i < n; ++i)
    {
        for (int j = i + 1; j <= n; ++j)
        {
            if (__gcd(a[i], a[j]) != 1)
            {
                dp[j] = (dp[j] + dp[i]) % mod;
            }
        }
    }
    cout << dp[n] << endl;
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