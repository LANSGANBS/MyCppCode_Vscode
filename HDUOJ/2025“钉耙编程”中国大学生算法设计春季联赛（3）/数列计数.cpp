#include <bits/stdc++.h>
// #include <bits/extc++.h>
using namespace std;
// using namespace __gnu_pbds;
#define endl '\n'
#define ture true
#define flase false
#define all(x) begin(x), end(x)
#define mem(a, x) memset(a, x, sizeof(a))
#define gcd(a, b) gcdint(a, b)
#define lcm(a, b) (a / gcd(a, b) * b)
#define sz(x) (int)x.size()
#define lowbit(x) (x & -x)
#define time(a, b) (abs((b - a) / CLOCKS_PER_SEC))
#define pb push_back
#define EPS 1e-7
using i64 = long long;
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

constexpr int inf = 0x7fffffff;
constexpr int N = 1.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

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
using Z = Zmod<MOD[0]>;

Z power(int n) { return mypow(Z(2), n); }

Z cnt(int a, int L) {
  if (L >= a) {
    int ones = __builtin_popcountll(a);
    return power(ones);
  }
  const int MMAX = 30;
  V<V<Z>> dp(MMAX + 1, vector<Z>(2, 0));
  dp[MMAX][1] = 1;
  for (int i = MMAX - 1; i >= 0; i--) {
    int bit_a = (a >> i) & 1;
    int bit_L = (L >> i) & 1;
    for (int j = 0; j <= 1; j++) {
      if (dp[i + 1][j] == 0) {
        continue;
      }
      dp[i][j & (bit_L == 0)] += dp[i + 1][j];
      if (bit_a == 1 && (!j || bit_L == 1)) {
        dp[i][j & (bit_L == 1)] += dp[i + 1][j];
      }
    }
  }
  return dp[0][0] + dp[0][1];
}

void solve() {
  int n;
  cin >> n;
  V<int> a(n), L(n);
  for (int i = 0; i < n; i++) {
    cin >> a[i];
  }
  for (int i = 0; i < n; i++) {
    cin >> L[i];
  }
  Z ans = 1;
  for (int i = 0; i < n; i++) {
    ans = (ans * cnt(a[i], L[i]));
  }
  cout << ans << endl;
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