#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

const int N = 2e5 + 7;
vector<int> e[N];
int dep[N << 1], st[N << 1 + 2][31], max_cnt[N << 1], d[N << 1], fa[N << 1], lg[N << 1];

void init(int n)
{
    for (int i = 1; i <= n; i++)
    {
        lg[i] = lg[i - 1] + (1 << lg[i - 1] == i);
    }

    for (int i = 0; i <= n; i++)
    {
        max_cnt[i] = fa[i] = dep[i] = d[i] = 0;
        e[i].clear();
        for (int j = 0; j <= 30; j++)
        {
            st[i][j] = 0;
        }
    }
}

int LCA(int x, int y)
{
    if (dep[x] < dep[y])
    {
        swap(x, y);
    }
    int d = dep[x] - dep[y];
    for (int i = 30; i >= 0; i--)
    {
        if (d >= (1 << i))
        {
            x = st[x][i];
            d -= (1 << i);
        }
    }
    if (x == y)
    {
        return x;
    }
    for (int i = 30; i >= 0; i--)
    {
        if (st[x][i] != st[y][i])
        {
            x = st[x][i];
            y = st[y][i];
        }
    }
    return st[x][0];
}

void dfs1(int x, int fx)
{
    dep[x] = dep[fx] + 1;
    st[x][0] = fx;
    for (int i = 1; i <= 30; i++)
    {
        st[x][i] = st[st[x][i - 1]][i - 1];
    }
    for (int i = 0; i < (int)e[x].size(); i++)
    {
        if (e[x][i] == fx)
        {
            continue;
        }
        dfs1(e[x][i], x);
    }
}

void dfs2(int u, int fx)
{
    fa[u] = fx;
    d[u] = d[fa[u]] + max_cnt[u];
    for (int i = 0; i < (int)e[u].size(); ++i)
    {
        if (e[u][i] == fx)
        {
            continue;
        }
        dfs2(e[u][i], u);
    }
}

void solve()
{
    int n;
    cin >> n;
    init(n);
    for (int i = 1; i <= n; i++)
    {
        cin >> max_cnt[i];
    }
    for (int i = 1; i < n; i++)
    {
        int u, v;
        cin >> u >> v;
        e[u].push_back(v);
        e[v].push_back(u);
    }
    dfs1(1, 0);
    dfs2(1, 0);
    int x, y;
    cin >> x >> y;
    int lc = LCA(x, y);
    int dis1 = d[x] - d[lc] + d[fa[y]] - d[lc] + max_cnt[lc];
    int dis2 = dep[x] - dep[lc] + dep[fa[y]] - dep[lc] + 1;
    if (dis1 == dis2 * 2)
    {
        cout << "yes" << endl;
    }
    else
    {
        cout << "no" << endl;
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