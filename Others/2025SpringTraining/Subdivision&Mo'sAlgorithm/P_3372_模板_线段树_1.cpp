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

int a[N];    // 存储原始数组
int st[N];   // 记录每个块的起始位置
int ed[N];   // 记录每个块的结束位置
int sum[N];  // 记录每个块的元素总和
int add[N];  // 记录每个块的懒标记（区间加法的增量）
int len;     // 每个块的大小（sqrt(n)）
int id[N];   // 记录每个元素所属的块号

// 处理区间 [l, r] 内所有元素加 k
void change(int l, int r, int k) {
  if (id[l] == id[r]) {  // 若 l 和 r 在同一个块
    for (int i = l; i <= r; i++) {
      a[i] += k;        // 更新原始数组
      sum[id[i]] += k;  // 更新该块的总和
    }
  } else {  // 若 l 和 r 在不同的块
    // 先处理 l 所在的块（部分覆盖）
    for (int i = l; i <= ed[l]; i++) {
      a[i] += k;
      sum[id[i]] += k;
    }
    // 处理 r 所在的块（部分覆盖）
    for (int i = st[r]; i <= r; i++) {
      a[i] += k;
      sum[id[i]] += k;
    }
    // 处理完全被 [l, r] 包含的完整块
    for (int i = id[l] + 1; i < id[r]; i++) {
      add[i] += k;  // 只修改懒标记
    }
  }
}

// 查询区间 [l, r] 内所有元素的和
int query(int l, int r) {
  int ans = 0;
  if (id[l] == id[r]) {  // 若 l 和 r 在同一个块
    for (int i = l; i <= r; i++) {
      ans += a[i] + add[id[i]];  // 需要加上懒标记的值
    }
  } else {  // 若 l 和 r 在不同的块
    // 处理 l 所在的块（部分覆盖）
    for (int i = l; i <= ed[l]; i++) {
      ans += a[i] + add[id[i]];
    }
    // 处理 r 所在的块（部分覆盖）
    for (int i = st[r]; i <= r; i++) {
      ans += a[i] + add[id[i]];
    }
    // 处理完全被 [l, r] 包含的完整块
    for (int i = id[l] + 1; i < id[r]; i++) {
      ans += sum[i] + add[i] * (ed[i] - st[i] + 1);  // 计算块和
    }
  }
  return ans;
}

void solve() {
  int n, m;
  cin >> n >> m;
  len = sqrt(n);  // 计算分块大小

  // 预处理分块信息
  for (int i = 1; i <= n; i++) {
    cin >> a[i];
    id[i] = (i - 1) / len + 1;      // 计算当前元素属于哪个块
    st[i] = (id[i] - 1) * len + 1;  // 计算块的起始位置
    ed[i] = min(id[i] * len, n);    // 计算块的结束位置
    sum[id[i]] += a[i];             // 计算块内元素的总和
  }

  while (m--) {
    int op, x, y, k;
    cin >> op >> x >> y;
    if (op == 1) {
      cin >> k;
      change(x, y, k);
    } else {
      cout << query(x, y) << endl;
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