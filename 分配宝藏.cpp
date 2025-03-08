#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N = 2.01e6;
const int M = 1.01e3;
const int mod = 1e9 + 7;
#define endl '\n'

void solve() {
  int n;
  cin >> n;
  int k = n / 2;
  int ans = (k % mod) * ((k + 1) % mod) % mod;
  cout << ans << endl;
}

signed main() {
  ios::sync_with_stdio(false);
  cin.tie(0), cout.tie(0);
  int T = 1;
  cin >> T;
  while (T--) solve();
}