#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e5;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

void solve()
{
    int w, b;
    cin >> w >> b;

    auto check = [&](int k)
    {
        int white = w, black = b;
        for (int i = k; i >= 1; --i)
        {
            if (white >= black)
            {
                if (white >= i)
                    white -= i;
                else if (black >= i)
                    black -= i;
                else
                    return false;
            }
            else
            {
                if (black >= i)
                    black -= i;
                else if (white >= i)
                    white -= i;
                else
                    return false;
            }
        }
        return true;
    };

    int left = 1, right = 2e9, ans = 0;
    while (left <= right)
    {
        int mid = left + (right - left) / 2;
        if (check(mid))
        {
            ans = mid;
            left = mid + 1;
        }
        else
        {
            right = mid - 1;
        }
    }
    cout << ans << endl;
}

signed main()
{
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    int T = 1;
    cin >> T;
    while (T--)
        solve();
}