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
#define gcd(a, b) gcdint(a, b)
#define lcm(a, b) (a / gcd(a, b) * b)
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
  constexpr Zmod operator-() const {
    Zmod val = norm(MOD - x);
    return val;
  }
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
  constexpr Zmod &operator%=(const int &i) { return x %= i, *this; }
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
using Z = Zmod<MOD[1]>;

Z power(int n) { return mypow(Z(2), n); }

// eg : Comb a;
// Z ans = a.C(100, 7);
// cout << ans << endl;

struct Comb {
  int n;
  vector<Z> _fac, _inv;

  Comb() : _fac{1}, _inv{0} {}
  Comb(int n) : Comb() { init(n); }
  void init(int m) {
    if (m <= n) return;
    _fac.resize(m + 1);
    _inv.resize(m + 1);
    for (int i = n + 1; i <= m; i++) {
      _fac[i] = _fac[i - 1] * i;
    }
    _inv[m] = _fac[m].inv();
    for (int i = m; i > n; i--) {
      _inv[i - 1] = _inv[i] * i;
    }
    n = m;
  }
  Z fac(int x) {
    if (x > n) init(x);
    return _fac[x];
  }
  Z inv(int x) {
    if (x > n) init(x);
    return _inv[x];
  }
  Z C(int x, int y) {
    if (x < 0 || y < 0 || x < y) return 0;
    return fac(x) * inv(y) * inv(x - y);
  }
  Z P(int x, int y) {
    if (x < 0 || y < 0 || x < y) return 0;
    return fac(x) * inv(x - y);
  }
} comb(1 << 21);

void solve() {
  Z n;
  cin >> n;
  V<Z> p(n.val());
  Z sum = 0, sumSq = 0;
  for (int i = 0; i < n.val(); i++) {
    cin >> p[i];
    sum += p[i];
    sumSq += p[i] * p[i];
  }
  Z S = (sum * sum - sumSq) / Z(2);
  Z ways = comb.fac(n.val()) * comb.fac(n.val() - 1) / power(n.val() - 1);
  cout << S << ' ' << ways << endl;
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