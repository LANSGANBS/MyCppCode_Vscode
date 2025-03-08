#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e2;
const int mod = 998244353;
#define endl '\n'
const int INF = 1e18;

void solve() {
  int n, m;
  cin >> n >> m;
  vector<vector<int>> t(n + 1, vector<int>(m + 1));
  vector<vector<int>> d(n + 1, vector<int>(m + 1));
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      cin >> t[i][j];
    }
  }
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      cin >> d[i][j];
    }
  }
  vector<vector<array<int, 2>>> dp(n + 1,
                                   vector<array<int, 2>>(m + 1, {INF, INF}));
  dp[1][1][0] = t[1][1];
  dp[1][1][1] = t[1][1] + d[1][1];
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      if (dp[i][j][0] != INF) {
        if (j + 1 <= m) {
          dp[i][j + 1][0] = min(dp[i][j + 1][0], dp[i][j][0] + t[i][j + 1]);
          dp[i][j + 1][1] =
              min(dp[i][j + 1][1], dp[i][j][0] + t[i][j + 1] + d[i][j + 1]);
        }
      }
      if (dp[i][j][1] != INF) {
        if (i + 1 <= n) {
          dp[i + 1][j][1] = min(dp[i + 1][j][1], dp[i][j][1] + t[i + 1][j]);
          dp[i + 1][j][0] =
              min(dp[i + 1][j][0], dp[i][j][1] + t[i + 1][j] + d[i + 1][j]);
        }
      }
    }
  }
  int ans = min(dp[n][m][1], dp[n][m][0] + d[n][m]);
  cout << ans << endl;
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  cin >> T;
  while (T--) solve();
}