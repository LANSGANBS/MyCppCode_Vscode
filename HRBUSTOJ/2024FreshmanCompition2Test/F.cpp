#include <bits/stdc++.h>
using namespace std;
#define int long long

const int MOD = 1e9 + 7;

int fastPow(int base, int exp)
{
    int res = 1;
    base %= MOD;
    while (exp > 0)
    {
        if (exp & 1)
            res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp >>= 1;
    }
    return res;
}

int modInverse(int num)
{
    return fastPow(num, MOD - 2);
}

void solve()
{
    int n, k;
    cin >> n >> k;
    vector<int> teams(n);

    int nonZerocnt = 0;
    int totalval = 0;

    for (int i = 0; i < n; i++)
    {
        cin >> teams[i];
        if (teams[i] > 0)
        {
            nonZerocnt++;
            totalval += teams[i];
            totalval %= MOD;
        }
    }

    int num = 0;
    int invN = modInverse(n);

    for (int i = 0; i < k; i++)
    {
        int curExp = (totalval + nonZerocnt * i) % MOD;
        num = (num + curExp * invN % MOD) % MOD;
    }

    cout << num << "\n";
}

signed main()
{
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    int tt = 1;
    cin >> tt;
    while (tt--)
    {
        solve();
    }
    return 0;
}