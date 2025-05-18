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
#define tcTV tcT, class... Ts

void unsyncIO() { cin.tie(0)->sync_with_stdio(0); }
void setPrec() { cout << fixed << setprecision(15); }
void setIO() { unsyncIO(), setPrec(); }

tcT > T gcd(const T &a, const T &b) { return b ? gcd(b, a % b) : a; }
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

tcTV > bool ckmin(T &a, const T &b, const Ts &...args) {
  bool changed = false;
  if (b < a) {
    a = b;
    changed = true;
  }
  (void)std::initializer_list<int>{
      ((args < a ? (a = args, changed = true) : 0), 0)...};
  return changed;
}

tcTV > bool ckmax(T &a, const T &b, const Ts &...args) {
  bool changed = false;
  if (a < b) {
    a = b;
    changed = true;
  }
  (void)std::initializer_list<int>{
      ((a < args ? (a = args, changed = true) : 0), 0)...};
  return changed;
}

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

struct DSU {
  vector<int> fa, p, e, f;

  DSU(int n) {
    fa.resize(n + 1);
    iota(fa.begin(), fa.end(), 0);
    p.resize(n + 1, 1);
    e.resize(n + 1);
    f.resize(n + 1);
  }
  int get(int x) {
    while (x != fa[x]) {
      x = fa[x] = fa[fa[x]];
    }
    return x;
  }
  bool merge(int x, int y) {  // 设x是y的祖先
    if (x == y) f[get(x)] = 1;
    x = get(x), y = get(y);
    e[x]++;
    if (x == y) return false;
    if (x < y) swap(x, y);  // 将编号小的合并到大的上
    fa[y] = x;
    f[x] |= f[y], p[x] += p[y], e[x] += e[y];
    return true;
  }
  bool same(int x, int y) { return get(x) == get(y); }
  bool F(int x) {  // 判断连通块内是否存在自环
    return f[get(x)];
  }
  int size(int x) {  // 输出连通块中点的数量
    return p[get(x)];
  }
  int E(int x) {  // 输出连通块中边的数量
    return e[get(x)];
  }
};

void solve() {
  int n, m;
  cin >> n >> m;
  V<V<int>> a(n, V<int>(m, 0));
  for (int i = 0; i < n; i++) {
    string s;
    cin >> s;
    assert(sz(s) == m);
    for (int j = 0; j < m; j++) {
      a[i][j] = s[j] - '0';
    }
  }
  auto id = [&](auto x, auto y) -> auto { return x * m + y; };
  DSU dsu(n * m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (a[i][j] == 0) {
        continue;
      }
      if (i > 0 and a[i - 1][j] == 1) {
        dsu.merge(id(i, j), id(i - 1, j));
      }
      if (j > 0 and a[i][j - 1] == 1) {
        dsu.merge(id(i, j), id(i, j - 1));
      }
    }
  }
  string ans;
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      if (a[i][j] == 0) {
        continue;
      }
      int cur = id(i, j);
      if (dsu.get(cur) != cur) {
        continue;
      }
      if (dsu.p[cur] <= n * 2) {
        continue;
      }
      if (dsu.p[cur] >= n * 15) {
        ans.pb('0');
      } else {
        ans.pb('1');
      }
    }
  }

  cout << ans << endl;
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