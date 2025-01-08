/*
 * URL: https://codeforces.com/contest/1918/problem/C#
 * @author: LANSGANBS
 * @Date: Created - 2024-01-31 00:11:38
 * @LastEditTime: 2024-01-31 13:32:10
 * 写完这道就去原
 */
#include <bits/stdc++.h>
#define endl '\n'
#define buff ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr)
#define REP(i, a, b) for (int i = a; i <= b; i++)
#define _REP(i, b, a) for (int i = b; i >= a; i--)
#define rep(i, a, b) for (int i = a; i < b; i++)
#define _rep(i, b, a) for (int i = b; i > a; i--)
#define all(x) x.begin(), x.end()
#define mem(a, x) memset(a, x, sizeof(a))
#define gcd(a, b) __gcd(a, b)
#define lcm(a, b) (a / gcd(a, b) * b)
/*
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("O3")
*/
#ifdef ONLINE_JUDGE
constexpr int N = 1e5 + 7;
#else
constexpr int N = 1e3 + 7;
#endif
using namespace std;
#define ll long long
#define int ll
#define ld long double
const int mod = 1e9 + 7;

struct node
{
    int x, y;
} a[N];

void solve()
{
    int x, y, r;
    cin >> x >> y >> r;
    bool is_high = true;
    int now = 0;
    for (int i = 60; i >= 0; i--)
    {
        bool px = (x & (1ll << i));
        bool py = (y & (1ll << i));
        if (px == py)
        {
            continue;
        }
        if (is_high)
        {
            now = (1ll << i) * (px - py);
            is_high = false;
        }
        else
        {
            if (r >= (1ll << i) && (px - py) * now > 0)
            {
                r -= (1ll << i);
                now -= (px - py) * (1ll << i);
            }
            else
            {
                now += (px - py) * (1ll << i);
            }
        }
    }
    if (now > 0)
    {
        cout << now << endl;
    }
    else
    {
        cout << -now << endl;
    }
    // wcout <<L" "<<endl;
    // cout << R"( )" << endl;
}

signed main()
{
    // setlocale(LC_ALL, "");
    buff;
    int t = 1;
    cin >> t;
    while (t--)
    {
        solve();
    }
    return 0;
}