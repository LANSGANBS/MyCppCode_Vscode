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

inline int gcd(int a, int b) { return b ? gcd(b, a % b) : a; }
inline i128 gcd128(i128 a, i128 b) { return b ? gcd(b, a % b) : a; }
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

void solve() {
  int N, M, Sx, Sy;
  cin >> N >> M >> Sx >> Sy;
  unordered_map<int, vector<int>> row;
  unordered_map<int, vector<int>> col;
  row.reserve(N * 2);
  row.max_load_factor(0.7f);
  col.reserve(N * 2);
  col.max_load_factor(0.7f);
  VP<int, int> houses(N);
  for (auto &p : houses) {
    cin >> p.first >> p.second;
  }
  for (auto &p : houses) {
    row[p.second].pb(p.first);
    col[p.first].pb(p.second);
  }
  for (auto &r : row) {
    sort(r.second.begin(), r.second.end());
  }
  for (auto &c : col) {
    sort(c.second.begin(), c.second.end());
  }
  unordered_set<int> visited;
  visited.reserve(N * 2);
  visited.max_load_factor(0.7f);
  auto hashf = [&](int x, int y) {
    unsigned int xx = (unsigned int)(x + 1000000000LL);
    unsigned int yy = (unsigned int)(y + 1000000000LL);
    return ((xx << 20) ^ yy);
  };
  int x = Sx, y = Sy, ans = 0;
  while (M--) {
    char d;
    int c;
    cin >> d >> c;
    int nx = x, ny = y;
    if (d == 'L') {
      nx -= c;
    } else if (d == 'R') {
      nx += c;
    } else if (d == 'U') {
      ny += c;
    } else {
      ny -= c;
    }
    if (d == 'L' || d == 'R') {
      int low = min(x, nx), high = max(x, nx);
      auto &vec = row[y];
      auto it1 = lower_bound(vec.begin(), vec.end(), low);
      auto it2 = upper_bound(vec.begin(), vec.end(), high);
      for (; it1 != it2; it1++) {
        int xx = *it1;
        unsigned int hval = hashf(xx, y);
        if (!visited.count(hval)) {
          visited.insert(hval);
          ans++;
        }
      }
    } else {
      int low = min(y, ny), high = max(y, ny);
      auto &vec = col[x];
      auto it1 = lower_bound(vec.begin(), vec.end(), low);
      auto it2 = upper_bound(vec.begin(), vec.end(), high);
      for (; it1 != it2; ++it1) {
        int yy = *it1;
        unsigned int hval = hashf(x, yy);
        if (!visited.count(hval)) {
          visited.insert(hval);
          ans++;
        }
      }
    }
    x = nx;
    y = ny;
  }
  cout << x << " " << y << " " << ans << "\n";
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