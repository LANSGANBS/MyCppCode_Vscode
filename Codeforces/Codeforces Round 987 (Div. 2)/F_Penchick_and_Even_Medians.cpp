#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long
const int N = 100;

void solve()
{
    int n;
    cin >> n;
    if (n == 1)
    {
        cout << 1 << endl;
        return;
    }
    pair<int, int> pr[N << 1];
    cout << "? " << n << ' ';
    for (int i = 0; i < n; i++)
    {
        cout << i + 1 << endl;
    }
    cout.flush();
    int cnt = 0;
    cin >> pr[cnt].first >> pr[cnt].second;
    if (pr[cnt].first > pr[cnt].second)
    {
        swap(pr[cnt].first, pr[cnt].second);
    }
    cnt++;
    cout << "? " << n / 2 << ' ';
    for (int i = 0; i < n / 2; i++)
    {
        cout << i + 1 << ' ';
    }
    cin >> pr[cnt].first >> pr[cnt].second;
    cnt++;
    if (pr[cnt].first > pr[cnt].second)
    {
        swap(pr[cnt].first, pr[cnt].second);
    }
    if (pr[cnt - 1].first == pr[cnt].first and pr[cnt - 1].second == pr[cnt].second) // 双相同
    {
        if (n == cnt - 1)
        {
            cout << pr[cnt - 1].first << ' ' << pr[cnt].first;
        }
        cin >> pr[cnt].first >> pr[cnt].second;
        cnt++;
    }
    else if (pr[cnt - 1].first == pr[cnt].first and pr[cnt - 1].second == pr[cnt].second or pr[cnt - 1].first == pr[cnt].second or pr[cnt - 1].second == pr[cnt].first) // 单相同
    {
    }
    else // 无相同
    {
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