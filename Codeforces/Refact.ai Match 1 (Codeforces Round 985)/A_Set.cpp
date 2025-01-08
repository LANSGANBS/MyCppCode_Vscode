#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

void solve()
{
    int l, r, k;
    cin >> l >> r >> k;
    int L = r / k;
    cout << max(0, L - l + 1) << endl;
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