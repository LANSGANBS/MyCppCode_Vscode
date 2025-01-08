#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

vector<int> son, e[N];
bool vis[N];

void solve()
{
    int n;
    cin >> n;
    son.clear();
    for (int i = 1; i <= n; i++)
    {
        e[i].clear();
    }
    for (int i = 1; i < n; i++)
    {
        int x, y;
        cin >> x >> y;
        e[x].push_back(y);
        e[y].push_back(x);
    }
    function<void(int, int)> dfs = [&](auto x, auto fa) -> auto
    {
        bool ok = false;
        for (auto &it : e[x])
        {
            if (it == fa)
            {
                continue;
            }
            ok = true;
            dfs(it, x);
        }
        if (!ok)
        {
            son.push_back(x);
        }
        vis[x] = (ok ^ 1);
    };
    dfs(1, 0);
    string s;
    cin >> s;
    s = " " + s;
    int e = 0, f = 0;
    for (auto &x : son)
    {
        if (s[x] == '?')
        {
            e++;
        }
    }
    for (int i = 2; i <= n; i++)
    {
        if (!vis[i] && s[i] == '?')
        {
            f++;
        }
    }
    if (s[1] == '?')
    {
        int cnta = 0, cntb = 0;
        for (auto &x : son)
        {
            if (s[x] == '0')
            {
                cnta++;
            }
            if (s[x] == '1')
            {
                cntb++;
            }
        }
        if (cnta == cntb)
        {
            cout << ((int)son.size() + f % 2) / 2 << endl;
        }
        else
        {
            cout << max(cnta, cntb) + e / 2 << endl;
        }
    }
    else
    {
        int cnt = 0;
        for (auto &x : son)
        {
            if (s[x] != '?' && s[x] != s[1])
            {
                cnt++;
            }
        }
        cout << cnt + (e + 1) / 2 << endl;
    }
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