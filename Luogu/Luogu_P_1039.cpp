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

constexpr int mod = 1e9 + 7;
constexpr int inf = 0x7fffffff;
constexpr int N = 1.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

int n, m, p, f, flag, id, bad, asd, good, two;
string name[10005], q[10005], aaa, ans;
char gh;
string day[8] = {"",
                 "Today is Monday.",
                 "Today is Tuesday.",
                 "Today is Wednesday.",
                 "Today is Thursday.",
                 "Today is Friday.",
                 "Today is Saturday.",
                 "Today is Sunday."};
map<string, int> l_or_w;
map<string, int> man;

struct node {
  string post, people;
} a[10005];

string fun(string x, string guilty, int as) {
  int ty = 0;
  string y = "", xx = "";
  for (int i = 0; i < x.size(); i++) {
    if (ty == 1) y += x[i];
    if (x[i] == ' ') ty = 1;
    if (ty == 0) xx += x[i];
  }
  int fo = 0;
  for (int i = 1; i <= n; i++) {
    if (xx == name[i]) {
      fo = 1;
      break;
    }
  }
  if (fo == 0) return "abcd";
  if (y == "is not guilty.") {
    if (xx == guilty)
      return "bad";
    else
      return "good";
  } else if (y == "is guilty.") {
    if (xx != guilty) {
      two = 1;
      return "bad";
    } else
      return "good";
  }
  return "abcd";
}

void solve() {
  cin >> n >> m >> p;
  for (int i = 1; i <= n; i++) cin >> name[i];
  getline(cin, q[0]);
  for (int i = 1; i <= p; i++) {
    getline(cin, q[i]);
    if (q[i][q[i].size() - 1] == '\n' || q[i][q[i].size() - 1] == '\r' ||
        q[i][q[i].size() - 1] == ' ') {
      q[i].erase(q[i].size() - 1, 1);
    }
    string m = "", o = "";
    f = 0;
    for (int j = 0; j < q[i].size(); j++) {
      if (f == 1 && q[i][j - 1] != ':') o += q[i][j];
      if (q[i][j] == ':') f = 1;
      if (f == 0) m += q[i][j];
    }
    a[i].people = m;
    a[i].post = o;
  }

  for (int i = 7; i >= 1; i--) {
    for (int j = 1; j <= n; j++) {
      for (int k = 1; k <= n; k++) {
        l_or_w[name[k]] = 0;
      }
      bad = good = flag = 0;
      for (int k = 1; k <= p; k++) {
        if (l_or_w[a[k].people] == 1) {
          id = 0;
          for (int l = 1; l <= 7; l++) {
            if (a[k].post == day[l]) {
              id = l;
              break;
            }
          }
          if (id != 0 && id == i) {
            flag = 1;
            break;
          }
          if (a[k].post == "I am guilty." && name[j] == a[k].people) {
            flag = 1;
            break;
          }
          if (a[k].post == "I am not guilty." && a[k].people != name[j]) {
            flag = 1;
            break;
          }
          if (fun(a[k].post, name[j], 0) == "good") {
            flag = 1;
            break;
          }
          continue;
        }
        if (l_or_w[a[k].people] == 2) {
          id = 0;
          for (int l = 1; l <= 7; l++) {
            if (a[k].post == day[l]) {
              id = l;
              break;
            }
          }
          if (id != 0 && id != i) {
            flag = 1;
            break;
          }
          if (a[k].post == "I am guilty." && name[j] != a[k].people) {
            two = 1;
            flag = 1;
            break;
          }
          if (a[k].post == "I am not guilty.") {
            if (a[k].people == name[j]) {
              flag = 1;
              break;
            }
          }
          if (fun(a[k].post, name[j], 1) == "bad") {
            flag = 1;
            break;
          }
          continue;
        }
        if (l_or_w[a[k].people] == 0) {
          if (a[k].post == "I am guilty.") {
            if (a[k].people == name[j]) {
              l_or_w[a[k].people] = 2;
              good++;
            } else {
              l_or_w[a[k].people] = 1;
              bad++;
            }
            continue;
          } else if (a[k].post == "I am not guilty.") {
            if (a[k].people != name[j]) {
              l_or_w[a[k].people] = 2;
              good++;
            } else {
              l_or_w[a[k].people] = 1;
              bad++;
            }
            continue;
          } else {
            id = 0;
            for (int l = 1; l <= 7; l++) {
              if (day[l] == a[k].post) {
                id = l;
                break;
              }
            }
            if (id != 0) {
              if (id == i) {
                l_or_w[a[k].people] = 2;
                good++;
              } else {
                l_or_w[a[k].people] = 1;
                bad++;
              }
              continue;
            }
            if (fun(a[k].post, name[j], 0) == "good") {
              l_or_w[a[k].people] = 2;
              good++;
            }
            if (fun(a[k].post, name[j], 0) == "bad") {
              l_or_w[a[k].people] = 1;
              bad++;
            }
          }
        }
      }
      if (flag == 1 || bad > m || good > n - m) {
        continue;
      }
      man[name[j]] = 1;
      ans = name[j];
    }
  }
  for (int i = 1; i <= n; i++) {
    if (man[name[i]] == 1) asd++;
  }
  if (asd > 1) {
    cout << "Cannot Determine" << endl;
  }
  if (asd == 1) {
    cout << ans << endl;
  }
  if (asd == 0) {
    cout << "Impossible" << endl;
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