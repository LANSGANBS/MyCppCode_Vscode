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
    vector<vector<int>> e(n);
    for (int i = 1; i < n; i++)
    {
        int p;
        cin >> p;
        p--;
        e[p].push_back(i);
    }
    function<int(int x)> dfs = [&](int x)
    {
        map<int, int> mp;
        int cnt = 0;
        for (auto y : e[x])
        {
            mp[dfs(y)]++;
        }
        if (mp.empty())
            return 1;
        if (e[x].size() == 1)
            return mp.begin()->first + 1;
        bool flag = true;
        int mx = 0;
        while (!mp.empty())
        {
            auto [x, cnt] = *mp.begin();
            mp.erase(mp.begin());
            mx = x;
            if (cnt >= 2)
                mp[x + 1] += cnt / 2;
            cnt %= 2;
            if (cnt && !mp.empty())
                flag = false;
        }
        if (!flag)
            mx++;
        return mx;
    };
    cout << dfs(0) - 1 << endl;
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