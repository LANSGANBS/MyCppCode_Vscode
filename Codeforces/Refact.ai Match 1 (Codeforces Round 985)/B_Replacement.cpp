#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

void solve()
{
    int n;
    cin >> n;
    string s, r;
    cin >> s >> r;

    int cnt0 = 0, cnt1 = 0;
    for (int i = 0; i < n; i += 1)
    {
        cnt0 += s[i] == '0';
        cnt1 += s[i] == '1';
    }
    for (int i = 0; i < n - 1; i++)
    {
        if (cnt0 == 0 || cnt1 == 0)
        {
            cout << "NO" << endl;
            return;
        }
        if (r[i] == '0')
            cnt1--;
        if (r[i] == '1')
            cnt0--;
    }
    cout << "YES" << endl;
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