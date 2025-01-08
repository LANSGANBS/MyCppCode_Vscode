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
int n, ans[N];

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

int query(int l, int r)
{
    cout << "? " << l << " " << r << endl;
    cout << flush;
    int res;
    cin >> res;
    return res;
}

void solve()
{
    int n;
    cin >> n;
    vector<int> ans(n + 1, -1);
    int pos = 0;
    for (int i = 2; i <= n; i++)
    {
        if (query(i - 1, i))
        {
            pos = i;
            break;
        }
    }
    if (!pos)
    {
        cout << "! IMPOSSIBLE" << endl;
        cout << flush;
        return;
    }
    int cs = query(1, pos);
    for (int i = 1; i < pos - cs; i++)
    {
        ans[i] = 1;
    }
    for (int i = pos - 1; i >= pos - cs; i--)
    {
        ans[i] = 0;
    }
    ans[pos] = 1;
    for (int i = pos + 1; i <= n; i++)
    {
        int temp = query(1, i);
        if (temp > cs)
        {
            ans[i] = 1;
        }
        else
        {
            ans[i] = 0;
        }
        cs = temp;
    }
    cout << "! ";
    for (int i = 1; i <= n; i++)
    {
        cout << ans[i];
    }
    cout << endl;
    cout << flush;
}

signed main()
{
    // ios::sync_with_stdio(false);
    // cin.tie(0), cout.tie(0);
    int T = 1;
    cin >> T;
    while (T--)
        solve();
}