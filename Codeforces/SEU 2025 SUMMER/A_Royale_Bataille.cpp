// C++
#include <bits/stdc++.h>
using namespace std;
const int MOD = 998244353;

int n, m;
vector<string> mp;

// 检查一行是否合法
bool valid_row(const string& row, int l1, int r1, int l2, int r2) {
  for (int i = 0; i < m; ++i) {
    char c = row[i];
    if (l1 <= i && i <= r1) {
      if (c != '?' && c != '1') return false;
    } else if (l2 <= i && i <= r2) {
      if (c != '?' && c != '2') return false;
    } else {
      if (c != '?' && c != '0') return false;
    }
  }
  // 1段和2段不能相邻或重叠
  if (r1 + 1 >= l2) return false;
  return true;
}

// 检查两行之间是否合法
bool valid_between(int l1a, int r1a, int l2a, int r2a, int l1b, int r1b,
                   int l2b, int r2b) {
  // 1段必须包含上一行的1段，2段必须包含上一行的2段
  if (l1b < l1a || r1b > r1a) return false;
  if (l2b > l2a || r2b < r2a) return false;
  // 1段和2段不能相邻或重叠
  if (r1b + 1 >= l2b) return false;
  return true;
}

int solve() {
  // dp[i][l1][r1][l2][r2] 当前行1段[l1,r1]，2段[l2,r2]的方案数
  // 由于每行只有一段1和一段2，且0段在中间，可以只用(l1,r1,l2,r2)描述
  vector<tuple<int, int, int, int>> configs;
  // 预处理所有合法的分段
  for (int l1 = 0; l1 < m; ++l1)
    for (int r1 = l1 - 1; r1 < m; ++r1)  // r1<l1表示没有1段
      for (int l2 = 0; l2 < m; ++l2)
        for (int r2 = l2 - 1; r2 < m; ++r2)  // r2<l2表示没有2段
          if (r1 < l2 || r2 < l1)            // 1段和2段不能重叠
            configs.emplace_back(l1, r1, l2, r2);

  unordered_map<uint64_t, int> dp, ndp;
  // 第一行初始化
  for (auto [l1, r1, l2, r2] : configs) {
    if (valid_row(mp[0], l1, r1, l2, r2)) {
      uint64_t key = ((uint64_t)l1 << 48) | ((uint64_t)r1 << 36) |
                     ((uint64_t)l2 << 24) | ((uint64_t)r2 << 12);
      dp[key] = 1;
    }
  }
  for (int i = 1; i < n; ++i) {
    ndp.clear();
    for (auto [l1a, r1a, l2a, r2a] : configs) {
      uint64_t keya = ((uint64_t)l1a << 48) | ((uint64_t)r1a << 36) |
                      ((uint64_t)l2a << 24) | ((uint64_t)r2a << 12);
      if (!dp.count(keya)) continue;
      for (auto [l1b, r1b, l2b, r2b] : configs) {
        if (!valid_row(mp[i], l1b, r1b, l2b, r2b)) continue;
        if (!valid_between(l1a, r1a, l2a, r2a, l1b, r1b, l2b, r2b)) continue;
        uint64_t keyb = ((uint64_t)l1b << 48) | ((uint64_t)r1b << 36) |
                        ((uint64_t)l2b << 24) | ((uint64_t)r2b << 12);
        ndp[keyb] = (ndp[keyb] + dp[keya]) % MOD;
      }
    }
    swap(dp, ndp);
  }
  int ans = 0;
  for (auto& [k, v] : dp) ans = (ans + v) % MOD;
  return ans;
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  int t;
  cin >> t;
  while (t--) {
    cin >> n >> m;
    mp.resize(n);
    for (int i = 0; i < n; ++i) cin >> mp[i];
    cout << solve() << '\n';
  }
  return 0;
}