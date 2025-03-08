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

struct kkk {
  int l;   // 左端点
  int r;   // 右端点
  int t;   // 此询问前修改数量
  int id;  // 询问编号
} q[N];

// 修改操作结构体
struct ttt {
  int id;   // 修改位置
  int val;  // 修改值
} c[N];

int v[N];     // 当前每个位置上的颜色（初始值读入后存入 v[]）
int vis[N];   // 统计每个颜色在当前区间出现次数
int sum = 0;  // 当前区间内不同颜色的数量

// 宏定义 add/del （ x 为颜色值 ）
#define add(x)                \
  {                           \
    if (++vis[x] == 1) sum++; \
  }
#define del(x)                \
  {                           \
    if (--vis[x] == 0) sum--; \
  }

// 全局当前区间区间边界，用于判断修改操作是否影响当前区间
int currentL, currentR;

// change 函数：对第 x 次修改操作进行“执行或撤销”
// 注意：若修改位置在当前区间 [currentL, currentR]
// 内，先撤销原来值的贡献，再加上修改值的贡献，最后通过 swap 更新
void change(int x) {
  int pos = c[x].id;  // 修改的位置
  if (pos >= currentL && pos <= currentR) {
    del(v[pos]);
    add(c[x].val);
  }
  swap(c[x].val, v[pos]);
}

void solve() {
  int n, m;
  cin >> n >> m;
  for (int i = 1; i <= n; i++) {
    cin >> v[i];
  }
  int totQuery = 0, totMod = 0;
  // 按时间顺序读入操作
  for (int i = 1; i <= m; i++) {
    char op;
    cin >> op;
    if (op == 'Q') {
      int L, R;
      cin >> L >> R;
      totQuery++;
      q[totQuery].l = L;
      q[totQuery].r = R;
      q[totQuery].t = totMod;  // 此查询前已有 totMod 次修改
      q[totQuery].id = totQuery;
    } else {
      int P, C_new;
      cin >> P >> C_new;
      totMod++;
      c[totMod].id = P;
      c[totMod].val = C_new;
    }
  }
  int blockSize = pow(n, 2.0 / 3);
  if (blockSize == 0) blockSize = 1;
  sort(q + 1, q + totQuery + 1, [=](const kkk &a, const kkk &b) {
    int ablock = a.l / blockSize, bblock = b.l / blockSize;
    if (ablock != bblock) return ablock < bblock;
    int rblockA = a.r / blockSize, rblockB = b.r / blockSize;
    if (rblockA != rblockB) return a.r < b.r;
    return a.t < b.t;
  });
  currentL = 1, currentR = 0;
  int currentT = 0;
  sum = 0;
  memset(vis, 0, sizeof(vis));
  vector<int> res(totQuery + 1, 0);

  // 依次处理每个查询
  for (int i = 1; i <= totQuery; i++) {
    // 调整时间（处理修改操作）
    while (currentT < q[i].t) {
      change(++currentT);
    }
    while (currentT > q[i].t) {
      change(currentT--);
    }
    while (currentR < q[i].r) {
      add(v[++currentR]);
    }
    while (currentR > q[i].r) {
      del(v[currentR--]);
    }
    while (currentL < q[i].l) {
      del(v[currentL++]);
    }
    while (currentL > q[i].l) {
      add(v[--currentL]);
    }
    // 记录当前区间答案，即不同颜色的数量
    res[q[i].id] = sum;
  }
  for (int i = 1; i <= totQuery; i++) {
    cout << res[i] << endl;
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