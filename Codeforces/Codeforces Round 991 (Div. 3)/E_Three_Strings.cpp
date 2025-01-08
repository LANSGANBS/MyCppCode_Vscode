#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'
#define sz(x) (int)x.size()

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

void solve() {
  string a, b, c;
  cin >> a >> b >> c;
  int n = sz(a), m = sz(b);
  vector<vector<int> > dp(n + 7, vector<int>(m + 7, mod));
  dp[0][0] = 0;
  int cs = mod;
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= m; j++) {
      if (i + j >= sz(c)) {
        continue;
      }
      if (i < n) {
        cs = dp[i][j] + (a[i] != c[i + j]);
        dp[i + 1][j] = min(dp[i + 1][j], cs);
      }
      if (j < m) {
        cs = dp[i][j];
        if (b[j] != c[i + j]) {
          cs++;
        }
        dp[i][j + 1] = min(dp[i][j + 1], cs);
      }
    }
  }
  cout << dp[n][m] << endl;
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  cin >> T;
  while (T--) solve();
}
