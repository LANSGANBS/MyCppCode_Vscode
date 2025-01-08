#include <bits/stdc++.h>
// #include <bits/extc++.h>
using namespace std;
// using namespace __gnu_pbds;
#define i64 long long
#define i128 __int128
#define endl '\n'
#define buff ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr)
#define ture true
#define flase false
#define pow power
#define all(x) begin(x), end(x)
#define mem(a, x) memset(a, x, sizeof(a))
#define gcd(a, b) __gcd(a, b)
#define mylcm(a, b) (a / gcd(a, b) * b)
#define sz(x) (int)x.size()
#define lowbit(x) (x & -x)
#define pb push_back
#define EPS 1e-7
#define int ll
#define ll long long
#define fr first
#define sc second
#define vi vector<int>

i64 ceilDiv(i64 n, i64 m)  // u
{
  if (n >= 0) {
    return (n + m - 1) / m;
  } else {
    return n / m;
  }
}

i64 floorDiv(i64 n, i64 m)  // d
{
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

int power(int a, i64 b, int p) {
  int res = 1;
  for (; b; b /= 2, a = 1LL * a * a % p) {
    if (b % 2) {
      res = 1LL * res * a % p;
    }
  }
  return res;
}

/*int power(int a, i64 b)
{
    int res = 1;
    for (; b; b /= 2, a = 1LL * a * a)
    {
        if (b % 2)
        {
            res = 1LL * res * a;
        }
    }
    return res;
}/*

/*std::ostream &operator<<(std::ostream &os, i128 n)
{
    std::string s;
    while (n)
    {
        s += '0' + n % 10;
        n /= 10;
    }
    std::reverse(s.begin(), s.end());
    return os << s;
}*/

/*i128 gcd(i128 a, i128 b)
{
    return b ? gcd(b, a % b) : a;
}*/

template <class T, class U>  // h
T fstTrue(T lo, T hi, U f) {
  hi++;
  assert(lo <= hi);
  while (lo < hi) {
    T mid = lo + (hi - lo) / 2;
    if (f(mid)) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  return lo;
}

template <class T, class U>  // q
T lstTrue(T lo, T hi, U f) {
  lo--;
  assert(lo <= hi);
  while (lo < hi) {
    T mid = lo + (hi - lo + 1) / 2;
    if (f(mid)) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }
  return lo;
}

const int mod = 1e9 + 7;
const int inf = 0x7fffffff;
constexpr int N = 2.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

void solve() {
  int n;
  cin >> n;
  vector<int> a(n);
  cin >> a;
  map<int, int> mp;
  for (auto &x : a) {
    mp[x]++;
  }
  debug(mp);
  int ans = mp.begin()->sc / 2;
  int c = 0;
  for (auto &[x, w] : mp) {
    c += w;
    ans = max(ans, c * (n - c));
  }
  cout << ans << endl;
}

signed main() {
  buff;
  int tt = 1;
  cin >> tt;
  while (tt--) {
    solve();
  }
  return 0;
}