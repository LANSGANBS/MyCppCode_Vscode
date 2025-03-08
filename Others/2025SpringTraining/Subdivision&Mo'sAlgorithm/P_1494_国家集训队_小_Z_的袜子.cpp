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

// a[]：存储输入数组，数据下标从 1 开始
int a[N];
// sum：全局变量，用于统计当前区间中“有效对数”（例如：当前数字的贡献总和）
int sum;
// ans1[]、ans2[]：用于存储每个查询答案的分子和分母，最终输出化简后的分数
int ans1[N], ans2[N];
// len：块长
int len;
// cnt[]：计数数组，统计某个数出现的次数（用于更新 sum 时的贡献）
int cnt[N];

// 结构体 query 用于存储每个查询区间的信息
struct query {
  int l, r, id;  // l：查询区间左端点，r：右端点，id：该查询在原输入中的编号
} q[N];

// add(x)：将数值 x 加入当前区间，并更新答案 sum
// 这里的更新逻辑为：每加入一个 x，就使得 sum 增加 cnt[x]（因为 x
// 与当前区间中所有等于 x 的数形成一对）
void add(int x) {
  sum += cnt[x];
  cnt[x]++;
}

// del(x)：将数值 x 从当前区间移除，并更新答案 sum
// 移除时，先减少 cnt[x]，然后将 sum 减去剩下的 cnt[x]（去掉 x
// 与剩余相同数字的配对数）
void del(int x) {
  cnt[x]--;
  sum -= cnt[x];
}

void solve() {
  int n, m;
  cin >> n >> m;
  len = n / sqrt(m);
  for (int i = 1; i <= n; i++) {
    cin >> a[i];
  }
  for (int i = 1; i <= m; i++) {
    cin >> q[i].l >> q[i].r;
    q[i].id = i;  // 记录原查询编号，便于后续答案按原序输出
  }
  // 先按照块编号排序，相同块内再按照右端点排序（奇数块按升序、偶数块按降序）
  sort(q + 1, q + m + 1, [](const query &lhs, const query &rhs) {
    if ((lhs.l - 1) / len != (rhs.l - 1) / len) return lhs.l < rhs.l;
    // 块编号奇偶不同的排序方式，能使区间扩展减少时间复杂度
    if (((lhs.l - 1) / len + 1) & 1) return lhs.r < rhs.r;
    return lhs.r > rhs.r;
  });
  // 然后依次将区间调整为每个查询的区间，同时更新答案
  for (int i = 1, l = 1, r = 0; i <= m; i++) {
    // 若查询区间为单个元素，则答案为 0/1（单个数无法形成配对）
    if (q[i].l == q[i].r) {
      ans1[q[i].id] = 0;
      ans2[q[i].id] = 1;
      continue;
    }
    // 扩大左边界：当当前 l 大于目标左端点时，依次将前面的元素加入区间
    while (l > q[i].l) {
      add(a[--l]);
    }
    // 扩大右边界：当当前 r 小于目标右端点时，依次将后面的元素加入区间
    while (r < q[i].r) {
      add(a[++r]);
    }
    // 缩小左边界：当当前 l 小于目标左端点时，依次将当前最左边的元素移除
    while (l < q[i].l) {
      del(a[l++]);
    }
    // 缩小右边界：当当前 r 大于目标右端点时，依次将当前最右边的元素移除
    while (r > q[i].r) {
      del(a[r--]);
    }
    // 处理当前区间答案
    // 若 sum 为 0，则区间内没有形成任何有效对数，答案记为 0/1
    if (sum == 0) {
      ans1[q[i].id] = 0;
      ans2[q[i].id] = 1;
      continue;
    }
    // 否则，将 sum 记为分子，区间内所有可能的配对数作为分母
    ans1[q[i].id] = sum;
    ans2[q[i].id] = (r - l + 1) * (r - l) / 2;
    // 化简分数：计算分子和分母的最大公约数，并分别除以该公约数
    int t = __gcd(ans1[q[i].id], ans2[q[i].id]);
    ans1[q[i].id] /= t;
    ans2[q[i].id] /= t;
  }
  for (int i = 1; i <= m; i++) {
    cout << ans1[i] << '/' << ans2[i] << endl;
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