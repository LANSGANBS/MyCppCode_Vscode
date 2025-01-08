#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

inline int read()
{
    int x = 0, y = 1;
    char ch = getchar();
    while (ch < '0' || ch > '9')
    {
        if (ch == '-')
        {
            y = -1;
        }
        ch = getchar();
    }
    while (ch >= '0' && ch <= '9')
    {
        x = x * 10 + (ch ^ 48);
        ch = getchar();
    }
    return x * y;
}

const int N = 400000;
int T, n, a[500010];
vector<int> d[500010], vec;

inline bool cmp(int a, int b)
{
    return d[a].size() < d[b].size();
}

void solve()
{
    n = read();
    for (int i = 1; i <= n; i++)
    {
        a[i] = read();
    }
    sort(a + 1, a + n + 1, cmp);
    if (d[a[2]].size() == 1)
    {
        cout << -1 << endl;
        return;
    }
    if (d[a[1]].size() == 1)
    {
        for (int j = 2; j <= n; j++)
        {
            if (a[j] < (a[1] << 1) || ((a[j] & 1) && ((a[j] - d[a[j]][0]) < (a[1] << 1))))
            {
                cout << -1 << endl;
                return;
            }
        }
        cout << a[1] << endl;
        return;
    }
    cout << 2 << endl;
    return;
}

signed main()
{
    cin.tie(0);
    cout.tie(0);
    ios::sync_with_stdio(0);
    for (int i = 2; i <= N; i++)
    {
        for (int j = i; j <= N; j += i)
        {
            d[j].push_back(i);
        }
        vec.push_back(i);
    }
    T = read();
    while (T--)
        solve();
}