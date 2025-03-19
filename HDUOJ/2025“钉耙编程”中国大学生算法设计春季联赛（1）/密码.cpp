#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'

void solve() {
  int n;
  cin >> n;
  map<int, int> mp;
  for (int i = 0; i < n; i++) {
    int u, v, w;
    cin >> u >> v >> w;
    if (!((w - v) % u)) {
      mp[(w - v) / u]++;
    }
    if (!((v - w) % u)) {
      mp[(v - w) / u]++;
    }
    if (!((w - u) % v)) {
      mp[(w - u) / v]++;
    }
    if (!((u - w) % v)) {
      mp[(u - w) / v]++;
    }
    if (!((v - u) % w)) {
      mp[(v - u) / w]++;
    }
    if (!((u - v) % w)) {
      mp[(u - v) / w]++;
    }
  }
  for (auto [x, y] : mp) {
    if (y % n == 0 and x >= 0) {
      cout << x << endl;
    }
  }
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  cin >> T;
  while (T--) solve();
}