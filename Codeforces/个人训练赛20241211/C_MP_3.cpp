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
  int n, I;
  cin >> n >> I;
  vector<int> a(n);
  for (int i = 0; i < n; ++i) {
    cin >> a[i];
  }
  sort(a.begin(), a.end());
  vector<int> val;
  vector<int> count;
  vector<int> pre(1, 0);
  int cnt = 1;
  for (int i = 1; i <= n; ++i) {
    if (i == n || a[i] != a[i - 1]) {
      val.push_back(a[i - 1]);
      count.push_back(cnt);
      pre.push_back(pre.back() + cnt);
      cnt = 1;
    } else {
      cnt++;
    }
  }
  int max_bits = (8 * I) / n;
  int K;
  if (max_bits >= 20) {
    K = val.size();
  } else {
    K = min((int)val.size(), 1LL << max_bits);
  }
  int cht = n;
  int tot_val = val.size();
  for (int lo = 0, hi = 0; hi < tot_val; hi++) {
    while (hi - lo + 1 > K) {
      lo++;
    }
    int unchanged = pre[hi + 1] - pre[lo];
    int changes = n - unchanged;
    cht = min(cht, changes);
  }
  cout << cht << endl;
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  // cin >> T;
  while (T--) solve();
}