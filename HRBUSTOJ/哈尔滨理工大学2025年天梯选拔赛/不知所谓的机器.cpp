#include <bits/stdc++.h>
// #include <bits/extc++.h>
using namespace std;
// using namespace __gnu_pbds;
#define endl '\n'
#define ture true
#define flase false
#define pow power
#define all(x) begin(x), end(x)
#define mem(a, x) memset(a, x, sizeof(a))
#define gcd(a, b) gcdint(a, b)
#define lcm(a, b) (a / gcd(a, b) * b)
#define sz(x) (int)x.size()
#define lowbit(x) (x & -x)
#define time(a, b) (abs((b - a) / CLOCKS_PER_SEC))
// double a = clock();
#define pb push_back
#define EPS 1e-7
#define i64 long long
#define i128 __int128
#define fr first
#define sc second
#define tcT template <class T
#define tcTU tcT, class U

void unsyncIO() { cin.tie(0)->sync_with_stdio(0); }
void setPrec() { cout << fixed << setprecision(15); }
void setIO() { unsyncIO(), setPrec(); }

inline int gcdint(int a, int b) { return b ? gcdint(b, a % b) : a; }
inline i128 gcd128(i128 a, i128 b) { return b ? gcd128(b, a % b) : a; }
inline int cdiv(int a, int b) { return a / b + ((a ^ b) > 0 && a % b); }
inline int fdiv(int a, int b) { return a / b - ((a ^ b) < 0 && a % b); }

tcT > using V = vector<T>;
tcTU > using PR = pair<T, U>;
tcTU > using MP = map<T, U>;
tcTU > using VP = vector<pair<T, U>>;
tcT > using pqg = priority_queue<T, vector<T>, greater<T>>;
tcT > using pql = priority_queue<T, vector<T>, less<T>>;

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

int n, m;
V<int> a(N + 1);

struct Node {
  int l, r, val;
  int lazyAnd, lazyOr, lazyXor;
} seg[4 * N];

void build(int k, int l, int r) {
  seg[k].l = l;
  seg[k].r = r;
  seg[k].lazyAnd = -1;
  seg[k].lazyOr = 0;
  seg[k].lazyXor = 0;
  if (l == r) {
    seg[k].val = a[l];
    return;
  }
  int mid = (l + r) / 2;
  build(2 * k, l, mid);
  build(2 * k + 1, mid + 1, r);
}

void apply(int k, int A, int B, int C) {
  seg[k].val = ((seg[k].val & A) | B) ^ C;
  seg[k].lazyAnd &= A;
  seg[k].lazyOr = (seg[k].lazyOr & A) | B;
  seg[k].lazyXor = (seg[k].lazyXor & A) ^ C;
}

void pushdown(int k) {
  if (seg[k].l == seg[k].r) {
    return;
  }
  int A = seg[k].lazyAnd, B = seg[k].lazyOr, C = seg[k].lazyXor;
  if (A == -1 && B == 0 && C == 0) {
    return;
  }
  apply(2 * k, A, B, C);
  apply(2 * k + 1, A, B, C);
  seg[k].lazyAnd = -1;
  seg[k].lazyOr = 0;
  seg[k].lazyXor = 0;
}

void update(int k, int L, int R, int A, int B, int C) {
  if (L <= seg[k].l && seg[k].r <= R) {
    apply(k, A, B, C);
    {
      return;
    }
  }
  pushdown(k);
  int mid = (seg[k].l + seg[k].r) / 2;
  debug(mid);
  if (L <= mid) {
    update(2 * k, L, R, A, B, C);
  }
  if (mid < R) {
    update(2 * k + 1, L, R, A, B, C);
  }
}

int query(int k, int p) {
  if (seg[k].l == seg[k].r) {
    return seg[k].val;
  }
  pushdown(k);
  int mid = (seg[k].l + seg[k].r) / 2;
  debug(mid);
  if (p <= mid) {
    return query(2 * k, p);
  } else {
    return query(2 * k + 1, p);
  }
}

void solve() {
  cin >> n >> m;
  for (int i = 1; i <= n; i++) {
    cin >> a[i];
  }
  V<int> ans;
  build(1, 1, n);
  while (m--) {
    int op;
    cin >> op;
    int l, r, x;
    if (op == 1) {
      cin >> l >> r >> x;
      update(1, l, r, x, 0, 0);
    } else if (op == 2) {
      cin >> l >> r >> x;
      update(1, l, r, ~x, x, 0);
    } else if (op == 3) {
      cin >> l >> r >> x;
      update(1, l, r, -1, 0, x);
    } else if (op == 4) {
      int p;
      cin >> p;
      ans.pb(query(1, p));
    }
  }
  for (int i = 0; i < sz(ans); i++) {
    if (i != sz(ans) - 1) {
      cout << ans[i] << endl;
    } else {
      cout << ans[i];
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