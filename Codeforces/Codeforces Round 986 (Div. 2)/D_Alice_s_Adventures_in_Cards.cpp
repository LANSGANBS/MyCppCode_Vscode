#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

const int N = 2e5 + 7;
#define all(x) begin(x), end(x)

void solve()
{
    int n;
    cin >> n;
    // qkj
    vector<vector<int>> a(3, vector<int>(n));
    vector<vector<int>> p(3, vector<int>(n));
    vector<int> f(3);
    vector<int> id(n);
    iota(all(id), 0);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> a[i][j];
        }
        sort(all(id), [&](int x, int y)
             { return a[i][x] < a[i][y]; });
        for (int j = 0; j < n; j++)
        {
            a[i][j] = id[j];
            p[i][a[i][j]] = j;
        }
    }
    vector<pair<int, int>> la(n, {-1, -1});
    queue<int> q;
    q.push(0);
    la[0] = {0, -1};
    while (!q.empty())
    {
        int x = q.front();
        q.pop();
        for (int j = 0; j < 3; j++)
        {
            int g = p[j][x];
            while (f[j] < g)
            {
                int nx = a[j][f[j]];
                if (la[nx].first == -1 && nx >= x)
                {
                    la[nx].first = x;
                    la[nx].second = j;
                    q.push(nx);
                }
                f[j]++;
            }
        }
    }
    if (la[n - 1].first == -1)
    {
        cout << "NO\n";
        return;
    }
    vector<pair<char, int>> ans;
    n--;
    while (n > 0)
    {
        char c;
        if (la[n].second == 0)
            c = 'q';
        else if (la[n].second == 1)
            c = 'k';
        else
            c = 'j';
        ans.push_back({c, n});
        n = la[n].first;
    }
    reverse(all(ans));
    cout << "YES\n";
    cout << ans.size() << endl;
    for (auto [c, p] : ans)
    {
        cout << c << " " << p + 1 << endl;
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