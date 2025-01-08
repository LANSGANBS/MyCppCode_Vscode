#include <bits/stdc++.h>
using namespace std;
#define LOCAL

void solve()
{
    int n, L, R;
    cin >> n >> L >> R;
    vector<int> errors(n);
    for (int i = 0; i < n; ++i)
    {
        cin >> errors[i];
    }

    vector<int> summ(n + 1, 0);
    for (int i = 0; i < n; ++i)
    {
        summ[i + 1] = summ[i] + errors[i];
    }

    int maxx = INT_MIN;
    int minn = INT_MAX;

    for (int len = L; len <= R; ++len)
    {
        for (int i = 0; i <= n - len; ++i)
        {
            int sum = summ[i + len] - summ[i];
            maxx = max(maxx, abs(sum));
            minn = min(minn, abs(sum));
        }
    }
    cout << maxx << " " << minn << endl;
}

int main()
{
#ifdef LOCAL
    freopen("C:\\MyCppCode_Vscode\\input.txt", "r", stdin);
#endif
    ios::sync_with_stdio(false);
    cin.tie(0);
    int T;
    cin >> T;
    while (T--)
    {
        solve();
    }
    return 0;
}