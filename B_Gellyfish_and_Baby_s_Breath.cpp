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

using i64 = long long;

template <class T>
constexpr T mypow(T n, i64 k) {
  T r = 1;
  for (; k; k /= 2, n *= n) {
    if (k % 2) {
      r *= n;
    }
  }
  return r;
}

template <int MOD>
struct Zmod {
  int x;
  Zmod(int x = 0) : x(norm(x % MOD)) {}
  Zmod(i64 x) : x(norm(x % MOD)) {}

  constexpr int norm(int x) const {
    if (x < 0) {
      x += MOD;
    }
    if (x >= MOD) {
      x -= MOD;
    }
    return x;
  }

  constexpr int val() const { return x; }

  constexpr Zmod operator-() const { return Zmod(norm(MOD - x)); }

  constexpr Zmod inv() const {
    assert(x != 0);
    return mypow(*this, MOD - 2);
  }

  friend constexpr auto &operator>>(istream &in, Zmod &j) {
    int v;
    in >> v;
    j = Zmod(v);
    return in;
  }

  friend constexpr auto &operator<<(ostream &o, const Zmod &j) {
    return o << j.val();
  }

  constexpr Zmod &operator++() {
    x = norm(x + 1);
    return *this;
  }

  constexpr Zmod &operator--() {
    x = norm(x - 1);
    return *this;
  }

  constexpr Zmod operator++(int) {
    Zmod tmp = *this;
    ++(*this);
    return tmp;
  }

  constexpr Zmod operator--(int) {
    Zmod tmp = *this;
    --(*this);
    return tmp;
  }

  constexpr Zmod &operator+=(const Zmod &i) {
    x = norm(x + i.x);
    return *this;
  }

  constexpr Zmod &operator-=(const Zmod &i) {
    x = norm(x - i.x);
    return *this;
  }

  constexpr Zmod &operator*=(const Zmod &i) {
    x = i64(x) * i.x % MOD;
    return *this;
  }

  constexpr Zmod &operator/=(const Zmod &i) { return *this *= i.inv(); }

  constexpr Zmod &operator%=(const int &i) {
    x %= i;
    return *this;
  }

  friend constexpr Zmod operator+(const Zmod i, const Zmod j) {
    return Zmod(i) += j;
  }

  friend constexpr Zmod operator-(const Zmod i, const Zmod j) {
    return Zmod(i) -= j;
  }

  friend constexpr Zmod operator*(const Zmod i, const Zmod j) {
    return Zmod(i) *= j;
  }

  friend constexpr Zmod operator/(const Zmod i, const Zmod j) {
    return Zmod(i) /= j;
  }

  friend constexpr Zmod operator%(const Zmod i, const int j) {
    return Zmod(i) %= j;
  }

  friend constexpr bool operator==(const Zmod i, const Zmod j) {
    return i.val() == j.val();
  }

  friend constexpr bool operator!=(const Zmod i, const Zmod j) {
    return i.val() != j.val();
  }

  friend constexpr bool operator<(const Zmod i, const Zmod j) {
    return i.val() < j.val();
  }

  friend constexpr bool operator>(const Zmod i, const Zmod j) {
    return i.val() > j.val();
  }
};

constexpr int MOD[] = {998244353, 1000000007};
using Z = Zmod<MOD[0]>;

Z power(int n) { return mypow(Z(2), n); }

V<Z> pow2 = {1};

void solve() {
  int n;
  cin >> n;
  V<int> p(n), q(n), pos1(n), pos2(n);
  for (int i = 0; i < n; i++) {
    cin >> p[i];
    pos1[p[i]] = i;
  }
  for (int i = 0; i < n; i++) {
    cin >> q[i];
    pos2[q[i]] = i;
  }

  if (sz(pow2) <= n) {
    int so = pow2.size();
    pow2.resize(n + 1);
    for (int i = so; i <= n; i++) {
      pow2[i] = pow2[i - 1] * 2;
    }
  }
  V<int> mx1(n), mx2(n);
  mx1[0] = p[0];
  mx2[0] = q[0];
  for (int i = 1; i < n; i++) {
    mx1[i] = max(mx1[i - 1], p[i]);
    mx2[i] = max(mx2[i - 1], q[i]);
  }

  V<Z> r(n);
  for (int i = 0; i < n; i++) {
    int e = max(mx1[i], mx2[i]);
    int s;
    if (mx1[i] > mx2[i]) {
      int j = pos1[e];
      int k = i - j;
      s = q[k];
    } else if (mx2[i] > mx1[i]) {
      int k = pos2[e];
      int j = i - k;
      s = p[j];
    } else {
      int j = pos1[e];
      int k = i - j;
      int max1 = (k >= 0 && k < n ? q[k] : -1);
      int x = pos2[e];
      int y = i - x;
      int max2 = (y >= 0 && y < n ? p[y] : -1);
      s = max(max1, max2);
    }
    Z val = pow2[e] + pow2[s];
    r[i] = val;
  }
  for (int i = 0; i < n; i++) {
    cout << r[i] << " ";
  }
  cout << endl;
}

signed main() {
  setIO();
  int tt = 1;
  cin >> tt;
  while (tt--) {
    solve();
  }
  return 0;
}