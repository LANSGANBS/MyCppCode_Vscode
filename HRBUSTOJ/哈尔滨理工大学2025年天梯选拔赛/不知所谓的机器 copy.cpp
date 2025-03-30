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

inline int applyFunc(int x, int A, int B, int C) { return ((x & A) | B) ^ C; }

void compose(int &lAnd, int &lOr, int &lXor, int A, int B, int C) {
  lAnd &= A;
  lOr = (lOr & A) | B;
  lXor = (lXor & A) ^ C;
}

struct bb {
  int l, r;
  V<int> a;
  int lazyAnd, lazyOr, lazyXor;
};

void solve() {
  int n, m;
  cin >> n >> m;
  V<int> a(n + 1);
  for (int i = 1; i <= n; i++) {
    cin >> a[i];
  }
  int B = max(1ll, (int)sqrt(n));
  int nb = (n + B - 1) / B;
  vector<bb> b(nb);
  for (int i = 0; i < nb; i++) {
    int L = i * B + 1, R = min(n, (i + 1) * B);
    b[i].l = L;
    b[i].r = R;
    b[i].a.resize(R - L + 1);
    for (int j = L; j <= R; j++) {
      b[i].a[j - L] = a[j];
    }
    b[i].lazyAnd = -1;
    b[i].lazyOr = 0;
    b[i].lazyXor = 0;
  }
  auto pushDown = [&](bb &a) {
    for (auto &x : a.a) {
      x = applyFunc(x, a.lazyAnd, a.lazyOr, a.lazyXor);
    }
    a.lazyAnd = -1;
    a.lazyOr = 0;
    a.lazyXor = 0;
  };
  auto updateRange = [&](int l, int r, int A, int B, int C) {
    for (int i = 0; i < nb; i++) {
      int L = b[i].l, R = b[i].r;
      if (r < L || R < l) {
        continue;
      }
      if (l <= L && R <= r) {
        compose(b[i].lazyAnd, b[i].lazyOr, b[i].lazyXor, A, B, C);
      } else {
        pushDown(b[i]);
        for (int j = 0; j < b[i].a.size(); j++) {
          int pos = b[i].l + j;
          if (l <= pos && pos <= r) {
            b[i].a[j] = applyFunc(b[i].a[j], A, B, C);
          }
        }
      }
    }
  };
  auto queryPos = [&](int p) -> int {
    for (int i = 0; i < nb; i++) {
      if (b[i].l <= p && p <= b[i].r) {
        int idx = p - b[i].l;
        return applyFunc(b[i].a[idx], b[i].lazyAnd, b[i].lazyOr, b[i].lazyXor);
      }
    }
    return -1;
  };
  V<int> ans;
  while (m--) {
    int op;
    cin >> op;
    int l, r, x;
    if (op == 1) {
      cin >> l >> r >> x;
      updateRange(l, r, x, 0, 0);
    } else if (op == 2) {
      cin >> l >> r >> x;
      updateRange(l, r, ~x, x, 0);
    } else if (op == 3) {
      cin >> l >> r >> x;
      updateRange(l, r, -1, 0, x);
    } else if (op == 4) {
      int p;
      cin >> p;
      ans.pb(queryPos(p));
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