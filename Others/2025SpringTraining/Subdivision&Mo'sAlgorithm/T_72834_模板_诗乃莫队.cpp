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
constexpr int M = 5.01e4;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

enum EventType { QUERY, MODIFY };

struct Query {
  int l, r, t, idx;
};

struct Modification {
  int pos;    // 修改位置
  int type;   // 1: 修改颜色, 2: 修改数字
  int prevC;  // 若是颜色修改，原来的颜色
  int nowC;
  int prevA;  // 若是数字修改，原来的数字
  int nowA;
};

int n, m, k, totQuery, totMod, curAns, ansArr[M], initA[M], curA[M], initC[M],
    curC[M];
V<int> compColors, colSumArr;
V<Query> queries;
V<Modification> mods;

inline int f(int x) {
  if (k == 1) {
    return x;
  } else if (k == 2) {
    return x * x;
  } else {
    return x * x * x;
  }
}

inline void removePos(int pos) {
  int color = curC[pos];
  int val = curA[pos];
  int oldVal = colSumArr[color];
  int newVal = oldVal - val;
  curAns -= (f(oldVal) - f(newVal));
  colSumArr[color] = newVal;
}

inline void addPos(int pos) {
  int color = curC[pos];
  int val = curA[pos];
  int oldVal = colSumArr[color];
  int newVal = oldVal + val;
  curAns += (f(newVal) - f(oldVal));
  colSumArr[color] = newVal;
}

void applyModification(int modIdx, int L, int R) {
  Modification &mod = mods[modIdx];
  int pos = mod.pos;
  if (L <= pos && pos <= R) {
    removePos(pos);
  }
  if (mod.type == 1) {
    curC[pos] = mod.nowC;
  } else {
    curA[pos] = mod.nowA;
  }
  if (L <= pos && pos <= R) {
    addPos(pos);
  }
}

void undoModification(int modIdx, int L, int R) {
  Modification &mod = mods[modIdx];
  int pos = mod.pos;
  if (L <= pos && pos <= R) {
    removePos(pos);
  }
  if (mod.type == 1) {
    curC[pos] = mod.prevC;
  } else {
    curA[pos] = mod.prevA;
  }
  if (L <= pos && pos <= R) {
    addPos(pos);
  }
}

void solve() {
  cin >> n >> m >> k;
  for (int i = 1; i <= n; i++) {
    cin >> initA[i];
    curA[i] = initA[i];
  }
  for (int i = 1; i <= n; i++) {
    cin >> initC[i];
    compColors.pb(initC[i]);
    curC[i] = initC[i];
  }
  int lastans = 0;
  for (int i = 0; i < m; i++) {
    char op;
    cin >> op;
    if (op == 'Q') {
      int l, r;
      cin >> l >> r;
      l ^= lastans;
      r ^= lastans;
      if (l > r) {
        swap(l, r);
      }
      queries.pb({l, r, totMod, totQuery});
      totQuery++;
    } else if (op == 'C') {
      int x, y;
      cin >> x >> y;
      x ^= lastans;
      y ^= lastans;
      compColors.pb(y);
      mods.pb({x, 1, curC[x], y, 0, 0});
      curC[x] = y;
      totMod++;
    } else {
      int x, y;
      cin >> x >> y;
      x ^= lastans;
      y ^= lastans;
      mods.pb({x, 2, 0, 0, curA[x], y});
      curA[x] = y;
      totMod++;
    }
  }
  sort(all(compColors));
  compColors.erase(unique(all(compColors)), compColors.end());
  int totColors = compColors.size();
  auto getColorId = [&](int c) -> int {
    return (int)(lower_bound(all(compColors), c) - compColors.begin()) + 1;
  };
  for (int i = 1; i <= n; i++) {
    initC[i] = getColorId(initC[i]);
    curC[i] = initC[i];
  }
  for (auto &mod : mods) {
    if (mod.type == 1) {
      mod.prevC = getColorId(mod.prevC);
      mod.nowC = getColorId(mod.nowC);
    }
  }
  for (int i = 1; i <= n; i++) {
    curA[i] = initA[i];
    curC[i] = initC[i];
  }
  colSumArr.assign(totColors + 1, 0);
  curAns = 0;
  int blockSize = max(1ll, pow(n, 2.0 / 3.0));
  sort(all(queries), [&](const Query &A, const Query &B) {
    int ablock = A.l / blockSize, bblock = B.l / blockSize;
    if (ablock != bblock) {
      return ablock < bblock;
    }
    int rblockA = A.r / blockSize, rblockB = B.r / blockSize;
    if (rblockA != rblockB) {
      return A.r < B.r;
    }
    return A.t < B.t;
  });
  for (int i = 1; i <= n; i++) {
    curA[i] = initA[i];
    curC[i] = initC[i];
  }
  int L = 1, R = 0, curT = 0;
  for (auto &q : queries) {
    while (curT < q.t) {
      applyModification(curT, L, R);
      curT++;
    }
    while (curT > q.t) {
      curT--;
      undoModification(curT, L, R);
    }
    while (R < q.r) {
      R++;
      addPos(R);
    }
    while (R > q.r) {
      removePos(R);
      R--;
    }
    while (L < q.l) {
      removePos(L);
      L++;
    }
    while (L > q.l) {
      L--;
      addPos(L);
    }
    ansArr[q.idx] = curAns;
    lastans = curAns;
  }
  for (int i = 0; i < totQuery; i++) {
    cout << ansArr[i] << endl;
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