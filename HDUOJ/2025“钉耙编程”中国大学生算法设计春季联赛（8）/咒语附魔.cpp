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

using ull = unsigned long long;
const int MOD1 = 1000000007;
const int MOD2 = 1000000009;
const int BASE = 91138233;

int add1(int a, int b) {
  a += b;
  if (a >= MOD1) a -= MOD1;
  return a;
}
int mul1(long long a, long long b) { return (a * b % MOD1); }
int add2(int a, int b) {
  a += b;
  if (a >= MOD2) a -= MOD2;
  return a;
}
int mul2(long long a, long long b) { return (a * b % MOD2); }

void solve() {
  int n, m;
  string A, B;
  cin >> n >> m >> A >> B;

  V<int> pw1(m + 1), pw2(m + 1);
  pw1[0] = pw2[0] = 1;
  for (int i = 1; i <= m; i++) {
    pw1[i] = mul1(pw1[i - 1], BASE);
    pw2[i] = mul2(pw2[i - 1], BASE);
  }

  V<int> h1(m + 1, 0), h2(m + 1, 0);
  for (int i = 0; i < m; i++) {
    int v = (B[i] - '0') + 1;
    h1[i + 1] = add1(mul1(h1[i], BASE), v);
    h2[i + 1] = add2(mul2(h2[i], BASE), v);
  }

  auto getHash = [&](int p, int len) {
    int x1 = h1[p + len] - mul1(h1[p], pw1[len]);
    if (x1 < 0) x1 += MOD1;
    int x2 = h2[p + len] - mul2(h2[p], pw2[len]);
    if (x2 < 0) x2 += MOD2;
    return PR<int, int>(x1, x2);
  };

  auto lcp = [&](int p, int q) {
    int lo = 0, hi = n;
    while (lo < hi) {
      int mid = midpoint(lo, hi);
      if (getHash(p, mid) == getHash(q, mid))
        lo = mid;
      else
        hi = mid - 1;
    }
    return lo;
  };

  auto better = [&](int p, int q) {
    int l = lcp(p, q);
    if (l >= n) {
      return false;
    }
    int cp = (A[l] - '0') ^ (B[p + l] - '0');
    int cq = (A[l] - '0') ^ (B[q + l] - '0');
    return cp > cq;
  };

  int best = 0;
  for (int p = 1; p <= m - n; p++) {
    if (better(p, best)) best = p;
  }

  int ans = 0;
  for (int i = 0; i < n; i++) {
    ans += ((A[i] - '0') ^ (B[best + i] - '0'));
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