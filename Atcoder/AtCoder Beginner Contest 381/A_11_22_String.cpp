/**
 *
 *
 *
 *                            iiiiii         iiiiiiiiii       iiiiiiiiiiiiii
 *                       iiiiiiiiiiiii   iiiiiii    iiii    iiiiiiiiiiiiiii                          ii   iiii
 *                    iiiiiiii     iiiiiiiii         iiii       iiii iii              iii          iiiiiiiiii
 *                 iiiiiii          iiiiii           iiii    iiii   ii           iiiiiiiiii      iiii iiii
 *               iiiiii            iiiii             iiii iiii        iii      iiii    iiiiiiiiiiiiiiiii  ii
 *             iiiiii            iiiiiii            iiiiiii       iiiiiiii   iii    iiiiiiiiiiiiii iii  iiii
 *           iiiiii             iiiiiii            iiiii   ii   iiii       iiiiiiiiiii iiii  iii iiii iiii      iii
 *          iiiii              iiiiiiii       ii        iiiii iiii    iiiiiiiii        iii iii  iii  iii  ii  iiii
 *        iiiiii              iiiiiiii      iiiii     iiiii iiiiiiiiiiiiiiii         iii  iii  ii  iii  iii iiii
 *       iiiii                 iiiiii     iiii     iiiiii iiiiiii    iii iii       iiii  ii   i   ii  iii  iii
 *     iiiiii                            iiii  iiiiiiiiiiiiiii       iii iiii   iiiii  iii  ii  iii  iii  ii
 *    iiiii                              iiiiiiii iiiiiiiiii       iiii   iiiiiiiii            ii  iii  ii
 *   iiiii                                     iiiiii  iiii      iiiii              iii      ii   ii  i
 * iiiiii                                  iiiiiiii   iiiii    iiiii                        ii  ii   ii
 * iiiii                                iiii  iiii    iiiiiiiiiiii                             ii
 *  iii                              iiii   iiii       iiiiiiii
 *                                iiiii   iiii
 *                              iiii     iiii
 *                            iiii    iiiii
 *                          iii     iiiii
 *                        iii     iiiii
 *                       iii   iiiiii
 *                       iiiiiiiii
 *                       iiiiii
 *
 *
 *
 **/
#include <bits/stdc++.h>
// #include <bits/extc++.h>
using namespace std;
// using namespace __gnu_pbds;
#define i64 long long
#define i128 __int128
#define endl '\n'
#define buff ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr)
#define debug cout << "----------------------------------------------" << endl
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
#define debug1(x)                         \
    {                                     \
        cerr << #x << " = " << x << "\n"; \
    };
#define debug2(x, y)                                                  \
    {                                                                 \
        cerr << #x << " = " << x << ", " << #y << " = " << y << "\n"; \
    };
#define debug3(x, y, z)                                                                           \
    {                                                                                             \
        cerr << #x << " = " << x << ", " << #y << " = " << y << ", " << #z << " = " << z << "\n"; \
    };
#define debug4(x, y, z, w)                                                                                                    \
    {                                                                                                                         \
        cerr << #x << " = " << x << ", " << #y << " = " << y << ", " << #z << " = " << z << ", " << #w << " = " << w << "\n"; \
    };

i64 ceilDiv(i64 n, i64 m) // u
{
    if (n >= 0)
    {
        return (n + m - 1) / m;
    }
    else
    {
        return n / m;
    }
}

i64 floorDiv(i64 n, i64 m) // d
{
    if (n >= 0)
    {
        return n / m;
    }
    else
    {
        return (n - m + 1) / m;
    }
}

template <typename T1, typename T2>
istream &operator>>(istream &in, pair<T1, T2> &a)
{
    return in >> a.first >> a.second;
}

template <typename T1>
istream &operator>>(istream &in, vector<T1> &a)
{
    for (auto &x : a)
    {
        in >> x;
    }
    return in;
}

template <typename T1, typename T2>
ostream &operator<<(ostream &out, const pair<T1, T2> &a)
{
    return out << a.first << ' ' << a.second;
}

template <typename T1, typename T2>
ostream &operator<<(ostream &out, const vector<pair<T1, T2>> &a)
{
    for (auto &x : a)
    {
        out << x << endl;
    }
    return out;
}

template <typename T1>
ostream &operator<<(ostream &out, const vector<T1> &a)
{
    int n = a.size();
    if (!n)
    {
        return out;
    }
    out << a[0];
    for (int i = 1; i < n; i++)
    {
        out << ' ' << a[i];
    }
    return out;
}

int power(int a, i64 b, int p)
{
    int res = 1;
    for (; b; b /= 2, a = 1LL * a * a % p)
    {
        if (b % 2)
        {
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

const int mod = 1e9 + 7;
const int inf = 0x7fffffff;
constexpr int N = 2.01e5;
constexpr int M = 2.01e3;

int n;
string s;

void solve()
{
    cin >> n;
    cin >> s;
    s = " " + s;
    if (n % 2 == 0)
    {
        cout << "No" << endl;
        return;
    }
    if (s[(n + 1) / 2] != '/')
    {
        cout << "No" << endl;
        return;
    }
    for (int i = 1; i < (n + 1) / 2; i++)
    {
        if (s[i] != '1')
        {
            cout << "No" << endl;
            return;
        }
    }
    for (int i = (n + 1) / 2 + 1; i <= n; i++)
    {
        if (s[i] != '2')
        {
            cout << "No" << endl;
            return;
        }
    }
    cout << "Yes";
}

signed main()
{
    buff;
    int tt = 1;
    // cin >> tt;
    while (tt--)
    {
        solve();
    }
    return 0;
}