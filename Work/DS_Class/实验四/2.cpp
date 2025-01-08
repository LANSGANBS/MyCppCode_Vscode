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

vector<int> g[N];
int n, m, u, v;
int vis[N];

void solve()
{
    memset(vis, 0, sizeof vis);
    cin >> n >> m;
    for (int i = 0; i < m; i++)
    {
        cin >> u >> v;
        g[u].push_back(v);
    }
    for (int i = 1; i <= n; i++)
    {
        sort(g[i].begin(), g[i].end());
    }
    function<void(int)> bfs = [&](int x) -> void
    {
        queue<int> q;
        q.push(x);
        vis[x] = 1;
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            cout << u << ' ';
            sort(g[u].begin(), g[u].end());
            for (auto c : g[u])
            {
                if (!vis[c])
                {
                    vis[c] = 1;
                    q.push(c);
                }
            }
        }
    };
    bfs(1);
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