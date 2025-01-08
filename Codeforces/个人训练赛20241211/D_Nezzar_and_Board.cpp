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

void solve() {
  int n, k;
  cin >> n >> k;
  vector<int> x(n);
  for (int i = 0; i < n; i++) {
    cin >> x[i];
  }
  int g = 0;
  for (int i = 1; i < n; i++) {
    g = std::__gcd(g, abs(x[i] - x[i - 1]));
  }
  int delta = abs(k - x[0]);
  if (delta % g == 0) {
    cout << "YES" << endl;
  } else {
    cout << "NO" << endl;
  }
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  cin >> T;
  while (T--) solve();
}