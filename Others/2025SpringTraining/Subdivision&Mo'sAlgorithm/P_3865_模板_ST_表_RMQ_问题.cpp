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
#define pb push_back
#define EPS 1e-7
// #define int ll
// #define ll long long
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

// st[]：每个块（分块）的起始下标
// ed[]：每个块的结束下标
// mx[]：每个块内的最大值（用于加速区间最大值查询）
// id[]：记录每个位置属于哪一个块（块编号）
// a[]：存放原始数组元素
// len：每个块的大小
int st[N], ed[N], mx[N], id[N], a[N], len;

int query() {
  int ans = 0;
  int l, r;
  cin >> l >> r;
  // 如果查询区间 l 与 r在同一块中，直接遍历区间求最大值
  if (id[l] == id[r]) {
    for (int i = l; i <= r; i++) {
      ans = max(ans, a[i]);
    }
  } else {
    // 分为三部分处理：
    // 1. l 所在块中，从下标 l 到该块的结束位置 ed[id[l]]
    for (int i = l; i <= ed[id[l]]; i++) {
      ans = max(ans, a[i]);
    }
    // 2. r 所在块中，从该块的起始位置 st[id[r]] 到 r
    for (int i = st[id[r]]; i <= r; i++) {
      ans = max(ans, a[i]);
    }
    // 3. 中间完全被 [l, r] 包含的块，直接取每个块的预处理最大值
    for (int i = id[l] + 1; i <= id[r] - 1; i++) {
      ans = max(ans, mx[i]);
    }
  }
  return ans;
}

void solve() {
  int n, m;
  cin >> n >> m;
  len = sqrt(n);
  for (int i = 1; i <= n; i++) {
    cin >> a[i];                // 读入第 i 个元素
    id[i] = (i - 1) / len + 1;  // 计算第 i 个元素所在的块编号
    st[id[i]] = (id[i] - 1) * len + 1;
    // 记录当前块的起始下标（对于同一块，每次赋值均相同）
    ed[id[i]] = id[i] * len;
    // 记录当前块的结束下标（可能超过 n，下文中后续处理时注意）
    // 更新当前块的最大值：若 i 是当前块的第一个元素，则直接赋值；否则更新最大值
    mx[id[i]] = (i == st[id[i]]) ? a[i] : max(mx[id[i]], a[i]);
  }
  // 注意：对于最后一块 ed[id] 可能超过 n，但本题数据保证 n 为正整数，
  // 实际查询时只遍历到 n 下标即可（若需要更严谨，可加 min(ed[id], n) 判断）
  // 处理 m 次查询操作，每次调用 query() 函数输出结果
  while (m--) {
    cout << query() << endl;
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