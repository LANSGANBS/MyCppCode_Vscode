#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 998244353;
const double eps = 1e-7;
#define endl '\n'

#ifdef LOCAL
#include "C:/Users/70510/Desktop/Others/algo/debug.h"
#else
#define debug(...) 42
#endif

int n;
double l;
double x[N], b[N];
double d[N];
int pre[N];

bool check(double mid) {
  for (int i = 0; i <= n; i++) {
    d[i] = 1e30;
  }
  d[0] = 0;
  pre[0] = -1;
  for (int i = 1; i <= n; i++) {
    for (int j = 0; j < i; j++) {
      double frustration = sqrt(abs(l - (x[i] - x[j])));
      double cost = d[j] + frustration - mid * b[i];
      if (d[i] > cost) {
        d[i] = cost;
        pre[i] = j;
      }
    }
  }
  return d[n] <= 0;
}

void solve() {
  cin >> n >> l;
  for (int i = 1; i <= n; i++) {
    cin >> x[i] >> b[i];
  }
  x[0] = 0;
  b[0] = 0;
  double lo = 0, hi = 1e6;
  while (hi - lo > eps) {
    double mid = (lo + hi) / 2;
    if (check(mid)) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  check(hi);
  vector<int> ans;
  for (int i = n; i != -1; i = pre[i]) {
    if (i != 0) {
      ans.push_back(i);
    }
  }
  reverse(ans.begin(), ans.end());
  for (int i = 0; i < (int)ans.size(); i++) {
    cout << ans[i] << " \n"[i == (int)ans.size() - 1];
  }
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  // cin >> T;
  while (T--) solve();
}