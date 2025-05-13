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
constexpr int mod = modulo[0];
constexpr int inf = 0x7fffffff;
constexpr int N = 1.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

int n;
V<int> adj[N], adj2[N];
PR<int, int> edges[N];
int fa1[N], sub[N], in[N], out[N], diff[N], cnt[N];
int timer1;
int subb[N], faa[N], index[N], order[N], ord, node, root;
bool used[N];
int id[N], col[N];
V<V<int>> g;

struct Item {
  int v, p, idx;
};

void solve() {
  cin >> n;
  for (int i = 1; i <= n; i++) {
    adj[i].clear();
  }
  for (int i = 1; i < n; i++) {
    int u, v;
    cin >> u >> v;
    edges[i] = {u, v};
    adj[u].pb(v);
    adj[v].pb(u);
  }

  timer1 = 0;
  for (int i = 1; i <= n + 1; i++) {
    diff[i] = 0;
  }
  vector<Item> st;
  st.reserve(n);
  st.pb({1, 0, 0});
  fa1[1] = 0;
  while (!st.empty()) {
    auto &it = st.back();
    int v = it.v, p = it.p;
    if (it.idx == 0) {
      fa1[v] = p;
      sub[v] = 1;
      in[v] = ++timer1;
    }
    if (it.idx < (int)adj[v].size()) {
      int u = adj[v][it.idx++];
      if (u == p) {
        continue;
      }
      st.pb({u, v, 0});
    } else {
      out[v] = timer1;
      if (p != 0) {
        int sz = sub[v];
        int oth = n - sz;
        if (sz <= oth) {
          diff[in[v]]++;
          diff[out[v] + 1]--;
        } else {
          diff[1]++;
          diff[n + 1]--;
          diff[in[v]]--;
          diff[out[v] + 1]++;
        }
        sub[p] += sub[v];
      }
      st.pop_back();
    }
  }
  for (int i = 1; i <= n; i++) {
    diff[i] += diff[i - 1];
  }
  for (int v = 1; v <= n; v++) {
    cnt[v] = diff[in[v]];
  }
  int bestl = inf, a = 1, b = 2;
  for (int i = 1; i < n; i++) {
    int x = edges[i].first;
    int y = edges[i].second;
    int a = min(x, y), b = max(x, y);
    int compb;
    if (fa1[b] == a) {
      compb = sub[b];
    } else {
      compb = n - sub[a];
    }
    int oth = n - compb;
    int orig_min = min(compb, oth);
    bool b_in_small = (compb <= oth);
    int l = orig_min + cnt[b] - (b_in_small ? 1 : 0);
    if (l < bestl) {
      bestl = l;
      a = a;
      b = b;
    }
  }
  node = n - 1;
  root = a;
  for (int i = 1; i <= n; i++) {
    adj2[i].clear();
  }
  for (int i = 1; i <= n; i++) {
    if (i == b) {
      continue;
    }
    for (auto &u : adj[i]) {
      if (u == b) {
        continue;
      }
      adj2[i].pb(u);
    }
  }
  for (auto &u : adj[b]) {
    if (u == a) {
      continue;
    }
    adj2[a].pb(u);
    adj2[u].pb(a);
  }
  ord = 0;
  for (int i = 1; i <= n; i++) {
    index[i] = 0;
    faa[i] = -1;
  }
  faa[root] = 0;
  V<int> stk;
  stk.reserve(node);
  stk.pb(root);
  while (!stk.empty()) {
    int v = stk.back();
    if (index[v] < (int)adj2[v].size()) {
      int u = adj2[v][index[v]++];
      if (u == faa[v]) {
        continue;
      }
      faa[u] = v;
      stk.pb(u);
    } else {
      stk.pop_back();
      order[ord++] = v;
    }
  }
  for (int i = 0; i < ord; i++) {
    int v = order[i];
    subb[v] = 1;
    for (int u : adj2[v]) {
      if (u == faa[v]) {
        continue;
      }
      subb[v] += subb[u];
    }
  }
  int cent = root, best = node;
  for (int idx = 0; idx < ord; idx++) {
    int v = order[idx];
    int h = node - subb[v];
    for (int u : adj2[v]) {
      if (u == faa[v]) {
        continue;
      }
      h = max(h, subb[u]);
    }
    if (h < best) {
      best = h;
      cent = v;
    }
  }

  for (int i = 1; i <= n; i++) {
    used[i] = false;
  }
  used[b] = true;
  used[cent] = true;
  g.clear();
  g.pb(V<int>(1, cent));
  int so = 1;
  for (int u : adj2[cent]) {
    if (used[u]) {
      continue;
    }
    g.pb(V<int>());
    stack<int> q;
    q.push(u);
    used[u] = true;
    while (!q.empty()) {
      int v = q.top();
      q.pop();
      g[so].pb(v);
      for (auto &w : adj2[v]) {
        if (!used[w]) {
          used[w] = true;
          q.push(w);
        }
      }
    }
    so++;
  }
  pqg<PR<int, int>> pq;
  for (int i = 0; i < sz(g); i++) {
    if (!g[i].empty()) pq.push({sz(g[i]), i});
  }
  for (int i = 1; i <= n; i++) {
    col[i] = 0;
  }
  col[b] = 0;
  int cur = 1;
  while (sz(pq) > 1) {
    auto p1 = pq.top();
    pq.pop();
    auto p2 = pq.top();
    pq.pop();
    int i1 = p1.sc, i2 = p2.sc;
    int v1 = g[i1].back();
    g[i1].pop_back();
    int v2 = g[i2].back();
    g[i2].pop_back();
    col[v1] = cur;
    col[v2] = cur;
    cur++;
    if (!g[i1].empty()) {
      pq.push({sz(g[i1]), i1});
    }
    if (!g[i2].empty()) {
      pq.push({sz(g[i2]), i2});
    }
  }
  int uout = a, vout = b;
  if (max(uout, vout) != b) {
    swap(uout, vout);
  }
  cout << uout << " " << vout << endl;
  for (int i = 1; i <= n; i++) {
    cout << col[i] << " \n"[i == n];
  }
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