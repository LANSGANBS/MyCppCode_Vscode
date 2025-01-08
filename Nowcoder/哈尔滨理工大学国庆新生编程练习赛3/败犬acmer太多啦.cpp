#include <bits/extc++.h>
#include <bits/stdc++.h>
using namespace std;
using namespace __gnu_pbds;
using i64 = long long;
#define endl '\n'
#define buff ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr)
#define debug cout << "----------------------------------------------" << endl
#define ture true
#define flase false
#define interesting int
#define all(x) begin(x), end(x)
#define mem(a, x) memset(a, x, sizeof(a))
#define gcd(a, b) __gcd(a, b)
#define lcm(a, b) (a / gcd(a, b) * b)
#define sz(x) (int)x.size()
#define lowbit(x) (x & -x)
#define pb push_back
#define EPS 1e-7
#define ll long long
#define int ll
#define ld long double
#define fr first
#define sc second
#define vi vector<int>
#define debug1(x) \
  { cerr << #x << " = " << x << "\n"; };
#define debug2(x, y) \
  { cerr << #x << " = " << x << ", " << #y << " = " << y << "\n"; };
#define debug3(x, y, z)                                                \
  {                                                                    \
    cerr << #x << " = " << x << ", " << #y << " = " << y << ", " << #z \
         << " = " << z << "\n";                                        \
  };
#define debug4(x, y, z, w)                                             \
  {                                                                    \
    cerr << #x << " = " << x << ", " << #y << " = " << y << ", " << #z \
         << " = " << z << ", " << #w << " = " << w << "\n";            \
  };

void yn(bool f) { cout << ((f == ture) ? "yes" : "no") << endl; }

void YN(bool F) { cout << ((F == ture) ? "YES" : "NO") << endl; }

i64 ceilDiv(i64 n, i64 m) {
  if (n >= 0) {
    return (n + m - 1) / m;
  } else {
    return n / m;
  }
}

i64 floorDiv(i64 n, i64 m) {
  if (n >= 0) {
    return n / m;
  } else {
    return (n - m + 1) / m;
  }
}

template <typename T1, typename T2>
istream &operator>>(istream &in, pair<T1, T2> &a) {
  return in >> a.first >> a.second;
}

template <typename T1>
istream &operator>>(istream &in, vector<T1> &a) {
  for (auto &x : a) {
    in >> x;
  }
  return in;
}

template <typename T1, typename T2>
ostream &operator<<(ostream &out, const pair<T1, T2> &a) {
  return out << a.first << ' ' << a.second;
}

template <typename T1, typename T2>
ostream &operator<<(ostream &out, const vector<pair<T1, T2>> &a) {
  for (auto &x : a) {
    out << x << endl;
  }
  return out;
}

template <typename T1>
ostream &operator<<(ostream &out, const vector<T1> &a) {
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

const int mod = 1e9 + 7;
constexpr int N = 2e5 + 7;
constexpr int M = 2e3 + 7;

void solve() {
  int m;
  cin >> m;
  vector<int> pw(1, 1), ans;
  for (int i = 1; i <= 11; i++) pw.push_back(pw[i - 1] * 3);
  while (m) {
    if (m == 1) {
      ans.push_back(0);
      break;
    }
    int k = 0;
    auto it = lower_bound(pw.begin(), pw.end(), m);
    if (*it > m) {
      it--;
    }
    m -= *it;
    ans.push_back(it - pw.begin());
  }
  cout << ans.size() << '\n';
  for (auto x : ans) {
    cout << x << ' ';
  }
}

signed main() {
  buff;
  int tt = 1;
  // cin >> tt;
  while (tt--) {
    solve();
  }
  return 0;
}