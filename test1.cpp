#include <bits/stdc++.h>
using namespace std;

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
  uniform_int_distribution<int> distA(1, 1000000000);

  int T = 10;
  int rem = 200;
  cout << T << "\n";
  for (int tc = 0; tc < T; tc++) {
    int max_n = min(rem - (T - tc - 1), 100000);
    int n = uniform_int_distribution<int>(1, max_n)(rng);
    rem -= n;

    int max_k = (n + 1) / 2;
    int k = uniform_int_distribution<int>(0, max_k)(rng);
    if (k == n) k = n - 1;

    cout << n << " " << k << "\n";
    for (int i = 0; i < n; i++) {
      cout << distA(rng) << (i + 1 == n ? "\n" : " ");
    }
  }
  return 0;
}