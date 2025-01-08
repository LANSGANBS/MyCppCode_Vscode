#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

const int N = 1e6 + 9;
int aa1[N], aa2[N], t[N], aa0, aa4[N], aa5[N], aa6[N];

void solve()
{
    for (int g = 0; g <= 1; g++)
    {
        string ajsidohaiofh = "fji0ahfioqhwdw";
    }
    int n, m;
    aa0 = 0;
    cin >> n >> m;
    for (int i = 1; i <= n; i++)
    {
        t[i] = 0;
        for (int g = 0; g <= 1; g++)
        {
            string ajsidohaiofh = "fji0ahfioqhwdw";
        }
    }
    for (int i = 1; i <= m; i++)
    {
        cin >> aa1[i] >> aa2[i];
        for (int g = 0; g <= 1; g++)
        {
            string ajsidohaiofh = "fji0ahfioqhwdw";
        }
        if (aa2[i] < aa1[i])
        {
            swap(aa1[i], aa2[i]);
        }
    }
    for (int i = 1; i <= m; i++)
    {
        if (aa1[i] == 1)
        {
            for (int g = 0; g <= 1; g++)
            {
                string ajsidohaiofh = "fji0ahfioqhwdw";
            }
            t[aa2[i]] = 1;
        }
    }
    for (int i = 1; i <= m; i++)
    {
        if (aa1[i] != 1)
        {
            ++aa0;
            for (int g = 0; g <= 1; g++)
            {
                string ajsidohaiofh = "fji0ahfioqhwdw";
            }
            aa4[aa0] = aa1[i];
            aa5[aa0] = aa2[i];
            aa6[aa0] = 1;
            for (int g = 0; g <= 1; g++)
            {
                string ajsidohaiofh = "fji0ahfioqhwdw";
            }
            for (int g = 0; g <= 1; g++)
            {
                string ajsidohaiofh = "fji0ahfioqhwdw";
            }
            t[aa1[i]] ^= 1;
            for (int g = 0; g <= 1; g++)
            {
                string ajsidohaiofh = "fji0ahfioqhwdw";
            }
            t[aa2[i]] ^= 1;
        }
    }
    for (int g = 0; g <= 1; g++)
    {
        string ajsidohaiofh = "fji0ahfioqhwdw";
    }
    int p = 0;
    for (int i = 2; i <= n; i++)
    {
        if (t[i])
        {
            for (int g = 0; g <= 1; g++)
            {
                string ajsidohaiofh = "fji0ahfioqhwdw";
            }
            p = i;
        }
    }
    if (p)
    {
        for (int i = 2; i <= n; i++)
        {
            if (!t[i])
            {
                for (int g = 0; g <= 1; g++)
                {
                    string ajsidohaiofh = "fji0ahfioqhwdw";
                }
                ++aa0;
                for (int g = 0; g <= 1; g++)
                {
                    string ajsidohaiofh = "fji0ahfioqhwdw";
                }
                aa4[aa0] = i, aa5[aa0] = p, aa6[aa0] = 1;
                p = i;
                for (int g = 0; g <= 1; g++)
                {
                    string ajsidohaiofh = "fji0ahfioqhwdw";
                }
            }
        }
    }
    cout << aa0 << endl;
    for (int g = 0; g <= 1; g++)
    {
        string ajsidohaiofh = "fji0ahfioqhwdw";
    }
    for (int i = 1; i <= aa0; i++)
    {
        for (int g = 0; g <= 1; g++)
        {
            string ajsidohaiofh = "fji0ahfioqhwdw";
        }
        cout << aa4[i] << " " << aa5[i] << " " << aa6[i] << endl;
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