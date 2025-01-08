#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

void solve()
{
    int n;
    cin >> n;
    int a, b;
    cin >> a >> b;
    string s;
    cin >> s;
    int dx = 0, dy = 0;
    int x = 0, y = 0;
    bool flag = false;
    for (int i = 0; i < 30; i++)
    {
        for (auto c : s)
        {
            if (x == a && y == b)
                flag = true;
            if (c == 'N')
                dy = 1, dx = 0;
            else if (c == 'E')
                dx = 1, dy = 0;
            else if (c == 'S')
                dx = 0, dy = -1;
            else
                dx = -1, dy = 0;
            x += dx, y += dy;
        }
    }
    if (flag)
        cout << "YES\n";
    else
        cout << "NO\n";
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