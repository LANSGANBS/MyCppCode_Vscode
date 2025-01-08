#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

#ifdef LOCAL
#include "C:/Users/70510/Desktop/Others/algo/debug.h"
#else
#define debug(...) 42
#endif

int x[N], y[N];

void solve() {
  int n;
  cin >> n;
  vector<double> xx;
  vector<double> yy;
  int totalPoints = n * 2;
  for (int i = 0; i < totalPoints; ++i) {
    int x, y;
    cin >> x >> y;
    if (x == 0) {
      yy.push_back(abs(y));
    } else if (y == 0) {
      xx.push_back(abs(x));
    }
  }
  sort(yy.begin(), yy.end());
  sort(xx.begin(), xx.end());
  double ans = 0.0;
  for (int i = 0; i < n; ++i) {
    double dis = sqrt(yy[i] * yy[i] + xx[i] * xx[i]);
    ans += dis;
  }
  cout << fixed << setprecision(15) << ans << endl;
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  cin >> T;
  while (T--) solve();
}