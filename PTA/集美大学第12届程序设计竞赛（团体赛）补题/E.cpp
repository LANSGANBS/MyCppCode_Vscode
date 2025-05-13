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
constexpr int MOD = modulo[1];
constexpr int inf = 0x7fffffff;
constexpr int N = 1.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

V<int> fac(int n) {
  V<int> fact(n + 1);
  fact[0] = 1;
  for (int i = 1; i <= n; i++) {
    fact[i] = fact[i - 1] * i;
  }
  return fact;
}

V<int> single(int a, const V<int> &fact, int n) {
  V<int> usage(n + 1, 0);
  for (int i = n; i >= 1; i--) {
    usage[i] = a / fact[i];
    a %= fact[i];
  }

  return usage;
}

V<int> cal(int l, int r, const V<int> &fact, int n) {
  V<int> total(n + 1, 0);
  for (int a = l; a <= r; a++) {
    V<int> usage = single(a, fact, n);
    for (int i = 1; i <= n; i++) {
      total[i] = (total[i] + usage[i]) % MOD;
    }
  }
  return total;
}

V<int> op(int l, int r, const V<int> &fact, int n) {
  V<int> res(n + 1, 0);
  for (int i = 1; i <= n; i++) {
    int cycle = fact[i];
    int com = (r / cycle) - ((l - 1) / cycle);
    int per = 0;
    for (int j = 0; j < cycle; j++) {
      V<int> usage = single(j, fact, n);
      per = (per + usage[i]) % MOD;
    }
    res[i] = (res[i] + (com * per) % MOD) % MOD;
    int l = (l - 1) / cycle * cycle + 1;
    for (int j = l; j < l; j++) {
      V<int> usage = single(j, fact, n);
      res[i] = (res[i] - usage[i] + MOD) % MOD;
    }
    int r = (r / cycle) * cycle;
    for (int j = r + 1; j <= r; j++) {
      V<int> usage = single(j, fact, n);
      res[i] = (res[i] + usage[i]) % MOD;
    }
  }
  return res;
}

void solve() {
  int n, q;
  cin >> n >> q;
  V<int> fact = fac(n);
  while (q--) {
    int l, r;
    cin >> l >> r;
    V<int> usage;
    if (r - l + 1 <= 1000000) {
      usage = cal(l, r, fact, n);
    } else {
      usage = op(l, r, fact, n);
    }
    for (int i = 1; i <= n; i++) {
      cout << usage[i] << ' ';
    }
    cout << endl;
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