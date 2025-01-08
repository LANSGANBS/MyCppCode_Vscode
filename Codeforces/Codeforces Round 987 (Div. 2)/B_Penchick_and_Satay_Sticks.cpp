#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

void solve()
{
    int n;
    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; i++)
        cin >> a[i], a[i]--;
    for (int i = 0; i < n; i++)
    {
        if (a[i] != i)
        {
            if (i + 1 >= n || abs(a[i] - a[i + 1]) != 1)
            {
                cout << "NO\n";
                return;
            }
            swap(a[i], a[i + 1]);
        }
        if (a[i] != i)
        {
            cout << "NO\n";
            return;
        }
    }
    cout << "YES\n";
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