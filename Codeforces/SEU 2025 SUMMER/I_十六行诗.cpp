#include <bits/extc++.h>
#include <bits/stdc++.h>
using namespace std;
using namespace __gnu_pbds;
#define endl '\n'
#define ture true
#define flase false
#define pow power
#define all(x) begin(x), end(x)
#define mem(a, x) memset(a, x, sizeof(a))
#define sz(x) (int)x.size()
#define lowbit(x) (x & -x)
#define time(a, b) (abs((b - a) / CLOCKS_PER_SEC))
// double a = clock();
#define pb push_back
#define EPS 1e-7
#define int ll
#define ll long long
#define i64 long long
#define i128 __int128
#define fr first
#define sc second
#define tcT template <class T
#define tcTU tcT, class U

void unsyncIO() { cin.tie(0)->sync_with_stdio(0); }
void setPrec() { cout << fixed << setprecision(15); }
void setIO() { unsyncIO(), setPrec(); }

tcT > T gcd(const T &a, const T &b) { return b ? gcd(b, a % b) : a; }
tcT > T lcm(const T &a, const T &b) { return a / gcd(a, b) * b; }
tcTU > T div(T a, T b, U flag) {
  if (flag) {
    return a / b + ((a ^ b) > 0 && a % b);
  } else {
    return a / b - ((a ^ b) < 0 && a % b);
  }
}

tcT > using V = vector<T>;
tcTU > using PR = pair<T, U>;
tcTU > using MP = map<T, U>;
tcTU > using VP = vector<pair<T, U>>;
tcT > using pql =
    __gnu_pbds::priority_queue<T, less<T>, __gnu_pbds::pairing_heap_tag>;
tcT > using pqg =
    __gnu_pbds::priority_queue<T, greater<T>, __gnu_pbds::pairing_heap_tag>;

tcTU > istream &operator>>(istream &in, pair<T, U> &a) {
  return in >> a.first >> a.second;
}

tcT > istream &operator>>(istream &in, vector<T> &a) {
  for (auto &x : a) {
    in >> x;
  }
  return in;
}

tcTU > ostream &operator<<(ostream &out, const pair<T, U> &a) {
  return out << a.first << ' ' << a.second;
}

tcTU > ostream &operator<<(ostream &out, const vector<pair<T, U>> &a) {
  for (auto &x : a) {
    out << x << endl;
  }
  return out;
}

tcT > ostream &operator<<(ostream &out, const vector<T> &a) {
  int n = a.size();
  if (!n) {
    return out;
  }
  out << a[0];
  for (int i = 1; i < n; i++) {
    out << ' ' << a[i];
  }
  return out;
}

std::ostream &operator<<(std::ostream &os, i128 n) {
  std::string s;
  while (n) {
    s += '0' + n % 10;
    n /= 10;
  }
  std::reverse(s.begin(), s.end());
  return os << s;
}

inline int power(int a, int b, int p = 1e9 + 7) {
  int res = 1;
  for (; b; b /= 2, a = 1LL * a * a % p) {
    if (b % 2) {
      res = 1LL * res * a % p;
    }
  }
  return res;
}

tcT > bool ckmin(T &a, const T &b) { return b < a ? a = b, 1 : 0; }
tcT > bool ckmax(T &a, const T &b) { return a < b ? a = b, 1 : 0; }

tcT > void remDup(vector<T> &v) {
  sort(all(v));
  v.erase(unique(all(v)), end(v));
}

tcTU > void erase(T &t, const U &u) {
  auto it = t.find(u);
  assert(it != end(t));
  t.erase(it);
}

tcTU > T fstTrue(T lo, T hi, U f) {
  hi++;
  assert(lo <= hi);
  while (lo < hi) {
    T mid = lo + (hi - lo) / 2;
    f(mid) ? hi = mid : lo = mid + 1;
  }
  return lo;
}

tcTU > T lstTrue(T lo, T hi, U f) {
  lo--;
  assert(lo <= hi);
  while (lo < hi) {
    T mid = lo + (hi - lo + 1) / 2;
    f(mid) ? lo = mid : hi = mid - 1;
  }
  return lo;
}

constexpr int modulo[] = {998244353, 1000000007};
constexpr int mod = modulo[0];
constexpr int inf = 0x7fffffff;
constexpr int N = 1.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

void solve() {
  int n;
  cin >> n;
  V<int> a(n + 1);
  for (int i = 1; i <= n; i++) {
    cin >> a[i];
  }
  V<int> tmp = a;
  remDup(tmp);
  int d = tmp.size();
  auto get = [&](int x) {
    return (lower_bound(tmp.begin() + 1, tmp.end(), x) - (tmp.begin()));
  };
  for (int i = 1; i <= n; i++) {
    a[i] = get(a[i]);
  }
  V<V<int>> nxt(n + 2, V<int>(d + 1, n + 1));
  for (int v = 1; v <= d; v++) {
    nxt[n + 1][v] = n + 1;
  }
  for (int i = n; i >= 1; i--) {
    for (int v = 1; v <= d; v++) {
      nxt[i][v] = nxt[i + 1][v];
    }
    nxt[i][a[i]] = i;
  }
  V<int> dp(n + 2, 0);
  for (int i = n; i >= 1; i--) {
    dp[i] = dp[i + 1];
    int j = nxt[i + 1][a[i]];
    if (j <= n) {
      for (int b = 1; b <= d; b++) {
        int k = nxt[j + 1][b];
        if (k > n) {
          continue;
        }
        int l = nxt[k + 1][b];
        if (l > n) {
          continue;
        }
        dp[i] = max(dp[i], 1 + dp[l + 1]);
      }
    }
    for (int j = i + 1; j <= n; j++) {
      int k = nxt[j + 1][a[i]];
      if (k > n) {
        continue;
      }
      int l = nxt[k + 1][a[j]];
      if (l > n) {
        continue;
      }
      dp[i] = max(dp[i], 1 + dp[l + 1]);
    }
    for (int j = i + 1; j <= n; j++) {
      int k = nxt[j + 1][a[j]];
      if (k > n) {
        continue;
      }
      int l = nxt[k + 1][a[i]];
      if (l > n) {
        continue;
      }
      dp[i] = max(dp[i], 1 + dp[l + 1]);
    }
  }
  cout << dp[1] << endl;
}

signed main() {
  setIO();
  int tt = 1;
  // cin >> tt;
  while (tt--) {
    solve();
  }
  return 0;
}