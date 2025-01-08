#include <bits/stdc++.h>
using namespace std;
using i64 = long long;
#define int long long
const int N = 1e5 + 10;
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

vector<pair<int, int>> ver[N];

void solve()
{
    int n, m;
    cin >> n >> m;
    for (int i = 0; i < m; i++)
    {
        int a, b, w;
        cin >> a >> b >> w;
        ver[a].push_back({b, w});
    }
    vector<int> dis(n + 1, 1E18);
    auto djikstra = [&](int s = 1) -> void
    {
        using PII = pair<int, int>;
        priority_queue<PII, vector<PII>, greater<PII>> q;
        q.emplace(0, s);
        dis[s] = 0;
        vector<int> vis(n + 1);
        while (!q.empty())
        {
            int x = q.top().second;
            q.pop();
            if (vis[x])
                continue;
            vis[x] = 1;
            for (auto [y, w] : ver[x])
            {
                if (dis[y] > dis[x] + w)
                {
                    dis[y] = dis[x] + w;
                    q.emplace(dis[y], y);
                }
            }
        }
    };
    djikstra(1);
    for (int i = 1; i <= n; i++)
    {
        cout << dis[i] << ' ';
    }
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