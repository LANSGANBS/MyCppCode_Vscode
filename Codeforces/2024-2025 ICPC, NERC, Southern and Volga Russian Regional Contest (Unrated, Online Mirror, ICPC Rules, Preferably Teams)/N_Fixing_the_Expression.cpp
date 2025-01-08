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

void solve()
{
    string s;
    cin >> s;
    if ((s[0] - '0') > (s[2] - '0') and s[1] != '>')
    {
        cout << s[0] << '>' << s[2] << endl;
    }
    else if ((s[0] - '0') < (s[2] - '0') and s[1] != '<')
    {
        cout << s[0] << '<' << s[2] << endl;
    }
    else if ((s[0] - '0') == (s[2] - '0') and s[1] != '=')
    {
        cout << s[0] << '=' << s[2] << endl;
    }
    else
    {
        cout << s << endl;
    }
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