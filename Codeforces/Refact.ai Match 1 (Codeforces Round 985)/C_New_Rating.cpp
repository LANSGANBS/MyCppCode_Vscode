#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define ll long long

void solve()
{
    int n;
    cin >> n;
    vector<int> a(n + 1);
    for (int i = 1; i <= n; i++)
        cin >> a[i];
    vector<vector<int>> dp(n + 3, vector<int>(3, 0));
    dp[0][1] = -1e9;
    dp[0][2] = -1e9;
    for (int i = 1; i <= n; i++)
    {
        int pr = 0;
        if (dp[i - 1][0] < a[i])
            pr = 1;

        /*
                if(dp[i-1][0] > a[i]) pr = -1;
                dp[i][0] = dp[i-1][0] + pr;
                dp[i][1] = max(dp[i-1][0] , dp[i-1][1]);
                int pr1 = 0 , pr2 = 0;
        //   if(dp[i-1][0] < a[i]) pr = 1;
        //         if(dp[i-1][0] > a[i]) pr = -1;
        //         dp[i][0] = dp[i-1][0] + pr;
        //         dp[i][1] = max(dp[i-1][0] , dp[i-1][1]);
        //         int pr1 = 0 , pr2 = 0;


        //         if(dp[i-1][1] < a[i]) pr1 = 1;
        //         if(dp[i-1][1] > a[i]) pr1 = -1;
        //         if(dp[i-1][2] < a[i]) pr2 = 1;

                if(dp[i-1][1] < a[i]) pr1 = 1;
                if(dp[i-1][1] > a[i]) pr1 = -1;

        //   if(dp[i-1][0] < a[i]) pr = 1;
        //         if(dp[i-1][0] > a[i]) pr = -1;
        //         dp[i][0] = dp[i-1][0] + pr;
        //         dp[i][1] = max(dp[i-1][0] , dp[i-1][1]);
        //         int pr1 = 0 , pr2 = 0;


        //         if(dp[i-1][1] < a[i]) pr1 = 1;
        //         if(dp[i-1][1] > a[i]) pr1 = -1;
        //         if(dp[i-1][2] < a[i]) pr2 = 1;

                if(dp[i-1][2] < a[i]) pr2 = 1;
                if(dp[i-1][2] > a[i]) pr2 = -1;
                dp[i][2] = max(dp[i-1][1] + pr1 ,dp[i-1][2] + pr2);
        */

        if (dp[i - 1][0] > a[i])
            pr = -1;
        dp[i][0] = dp[i - 1][0] + pr;
        dp[i][1] = max(dp[i - 1][0], dp[i - 1][1]);
        int pr1 = 0, pr2 = 0;
        //   if(dp[i-1][0] < a[i]) pr = 1;
        //         if(dp[i-1][0] > a[i]) pr = -1;
        //         dp[i][0] = dp[i-1][0] + pr;
        //         dp[i][1] = max(dp[i-1][0] , dp[i-1][1]);
        //         int pr1 = 0 , pr2 = 0;

        //         if(dp[i-1][1] < a[i]) pr1 = 1;
        //         if(dp[i-1][1] > a[i]) pr1 = -1;
        //         if(dp[i-1][2] < a[i]) pr2 = 1;

        if (dp[i - 1][1] < a[i])
            pr1 = 1;
        if (dp[i - 1][1] > a[i])
            pr1 = -1;

        //   if(dp[i-1][0] < a[i]) pr = 1;
        //         if(dp[i-1][0] > a[i]) pr = -1;
        //         dp[i][0] = dp[i-1][0] + pr;
        //         dp[i][1] = max(dp[i-1][0] , dp[i-1][1]);
        //         int pr1 = 0 , pr2 = 0;

        //         if(dp[i-1][1] < a[i]) pr1 = 1;
        //         if(dp[i-1][1] > a[i]) pr1 = -1;
        //         if(dp[i-1][2] < a[i]) pr2 = 1;

        if (dp[i - 1][2] < a[i])
            pr2 = 1;
        if (dp[i - 1][2] > a[i])
            pr2 = -1;
        dp[i][2] = max(dp[i - 1][1] + pr1, dp[i - 1][2] + pr2);
    }
    cout << max(dp[n][1], dp[n][2]) << endl;
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