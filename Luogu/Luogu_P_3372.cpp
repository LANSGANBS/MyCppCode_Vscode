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

inline int power(int a, i64 b, int p = 1e9 + 7) {
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

constexpr int mod = 1e9 + 7;
constexpr int inf = 0x7fffffff;
constexpr int N = 1.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

template <class T>
struct Segt_ {
  struct node {
    int l, r;
    T w, add, mul = 1;  // 注意初始赋值
  };
  vector<T> w;
  vector<node> t;

  Segt_(int n) {
    w.resize(n + 1);
    t.resize((n << 2) + 1);
    build(1, n);
  }
  Segt_(vector<int> in) {
    int n = in.size() - 1;
    w.resize(n + 1);
    for (int i = 1; i <= n; i++) {
      w[i] = in[i];
    }
    t.resize((n << 2) + 1);
    build(1, n);
  }
  void pushdown(node &p, T add, T mul) {  // 在此更新下递函数
    p.w = p.w * mul + (p.r - p.l + 1) * add;
    p.add = p.add * mul + add;
    p.mul *= mul;
  }
  void pushup(node &p, node &l, node &r) {  // 在此更新上传函数
    p.w = l.w + r.w;
  }
#define GL (k << 1)
#define GR (k << 1 | 1)
  void pushdown(int k) {  // 不需要动
    pushdown(t[GL], t[k].add, t[k].mul);
    pushdown(t[GR], t[k].add, t[k].mul);
    t[k].add = 0, t[k].mul = 1;
  }
  void pushup(int k) {  // 不需要动
    pushup(t[k], t[GL], t[GR]);
  }
  void build(int l, int r, int k = 1) {
    if (l == r) {
      t[k] = {l, r, w[l]};
      return;
    }
    t[k] = {l, r};
    int mid = (l + r) / 2;
    build(l, mid, GL);
    build(mid + 1, r, GR);
    pushup(k);
  }
  void modify(int l, int r, T val, int k = 1) {  // 区间修改
    if (l <= t[k].l && t[k].r <= r) {
      t[k].w += (t[k].r - t[k].l + 1) * val;
      t[k].add += val;
      return;
    }
    pushdown(k);
    int mid = (t[k].l + t[k].r) / 2;
    if (l <= mid) modify(l, r, val, GL);
    if (mid < r) modify(l, r, val, GR);
    pushup(k);
  }
  void modify2(int l, int r, T val, int k = 1) {  // 区间修改
    if (l <= t[k].l && t[k].r <= r) {
      t[k].w *= val;
      t[k].add *= val;
      t[k].mul *= val;
      return;
    }
    pushdown(k);
    int mid = (t[k].l + t[k].r) / 2;
    if (l <= mid) modify2(l, r, val, GL);
    if (mid < r) modify2(l, r, val, GR);
    pushup(k);
  }
  T ask(int l, int r, int k = 1) {  // 区间询问，不合并
    if (l <= t[k].l && t[k].r <= r) {
      return t[k].w;
    }
    pushdown(k);
    int mid = (t[k].l + t[k].r) / 2;
    T ans = 0;
    if (l <= mid) ans += ask(l, r, GL);
    if (mid < r) ans += ask(l, r, GR);
    return ans;
  }
#undef GL
#undef GR
};

void solve() {
  int n, m;
  cin >> n >> m;
  V<int> a(n + 1);
  for (int i = 1; i <= n; i++) {
    cin >> a[i];
  }
  Segt_<int> segt(a);
  int x, y, k;
  char ch;
  while (m--) {
    cin >> ch;
    if (ch == '1') {
      cin >> x >> y >> k;
      segt.modify(x, y, k);
    } else {
      cin >> x >> y;
      cout << segt.ask(x, y) << endl;
    }
  }
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