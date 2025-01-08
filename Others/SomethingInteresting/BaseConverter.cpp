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

const int mod = 1e9 + 7;
const int inf = 0x7fffffff;
constexpr int N = 2.01e6;
constexpr int M = 2.01e3;

#ifdef LOCAL
#include <C:/Users/70510/Desktop/Others/algo/debug.h>
#else
#define debug(...) 42
#endif

// BaseConverter::xtox(a, b, s)  // 将s的a进制数字转换为b进制数字

struct BaseConverter {
  static int charToValue(char ch) {
    if (isdigit(ch))
      return ch - '0';
    else if (isalpha(ch))
      return toupper(ch) - 'A' + 10;
    else
      return -1;
  }

  static char valueToChar(int val) {
    if (val >= 0 && val <= 9)
      return val + '0';
    else if (val >= 10 && val < 36)
      return val - 10 + 'A';
    else
      return '?';
  }
  static std::string xtox(int a, int b, const std::string &c) {
    long long value = 0;
    for (char ch : c) {
      int digit = charToValue(ch);
      if (digit < 0 || digit >= a) {
        throw std::invalid_argument("Error");
      }
      value = value * a + digit;
    }
    if (value == 0) return "0";
    std::string result;
    while (value > 0) {
      int remainder = value % b;
      result += valueToChar(remainder);
      value /= b;
    }
    std::reverse(result.begin(), result.end());
    return result;
  }
};

void solve() {
  cout << "+ - * /" << endl;
  string s;
  cin >> s;
  cout << "a进制的s1 b进制的s2 c进制的ans" << endl;
  int a, b, c;
  string s1, s2;
  cin >> a >> s1 >> b >> s2 >> c;
  if (s == "+") {
    cout << BaseConverter::xtox(
                10, c,
                (to_string(stoll(BaseConverter::xtox(a, 10, s1)) +
                           stoll(BaseConverter::xtox(b, 10, s2)))))
         << endl;
  } else if (s == "-") {
    cout << BaseConverter::xtox(
                10, c,
                (to_string(stoll(BaseConverter::xtox(a, 10, s1)) -
                           stoll(BaseConverter::xtox(b, 10, s2)))))
         << endl;
  } else if (s == "*") {
    cout << BaseConverter::xtox(
                10, c,
                (to_string(stoll(BaseConverter::xtox(a, 10, s1)) *
                           stoll(BaseConverter::xtox(b, 10, s2)))))
         << endl;
  } else if (s == "/") {
    cout << BaseConverter::xtox(
                10, c,
                (to_string(stoll(BaseConverter::xtox(a, 10, s1)) /
                           stoll(BaseConverter::xtox(b, 10, s2)))))
         << endl;
  }
}

signed main() {
  // setIO();
  int tt = 1;
  // cin >> tt;
  while (tt--) {
    solve();
  }
  return 0;
}