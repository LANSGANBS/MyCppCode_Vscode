#include <bits/stdc++.h>
#include <bits/extc++.h>
using namespace std;
using namespace __gnu_pbds;
using i64 = long long;
/*using i128 = __int128;*/
#define endl '\n'
#define buff ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr)
#define debug cout << "----------------------------------------------" << endl
#define ture true
#define flase false
#define pow power
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

int power(int a, i64 b)
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
}

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
constexpr int N = 2e5 + 7;
constexpr int M = 2e3 + 7;

template <int Mod>
class ModInt
{
private:
    int value;

public:
    static constexpr int mod = Mod;

    ModInt() : value(0) {}

    ModInt(int v)
    {
        value = v % Mod;
        if (value < 0)
        {
            value += Mod;
        }
    }

    int get_value() const
    {
        return value;
    }

    ModInt operator+(const ModInt &other) const
    {
        return ModInt(*this) += other;
    }

    ModInt &operator+=(const ModInt &other)
    {
        value = (value + other.value) % Mod;
        if (value < 0)
        {
            value += Mod;
        }
        return *this;
    }

    ModInt operator-(const ModInt &other) const
    {
        return ModInt(*this) -= other;
    }

    ModInt &operator-=(const ModInt &other)
    {
        value = (value - other.value + Mod) % Mod;
        return *this;
    }

    ModInt operator*(const ModInt &other) const
    {
        return ModInt(value * other.value) %= Mod;
    }

    ModInt &operator*=(const ModInt &other)
    {
        value = (value * other.value) % Mod;
        return *this;
    }

    ModInt operator/(const ModInt &other) const
    {
        return ModInt(*this) /= other;
    }

    ModInt &operator/=(const ModInt &other)
    {
        value = (value * other.inverse()) % Mod;
        return *this;
    }

    ModInt inverse() const
    {
        return power(*this, Mod - 2);
    }

    static ModInt power(ModInt base, int exponent)
    {
        ModInt result(1);
        ModInt baseValue = base;
        while (exponent > 0)
        {
            if (exponent % 2 == 1)
            {
                result *= baseValue;
            }
            baseValue *= baseValue;
            exponent /= 2;
        }
        return result;
    }

    friend std::ostream &operator<<(std::ostream &os, const ModInt &m)
    {
        os << m.value;
        return os;
    }

    ModInt &operator%=(int m)
    {
        value %= m;
        return *this;
    }

private:
    static int inverseMod(int a, int m)
    {
        a = a % m;
        for (int x = 1; x < m; x++)
        {
            if ((a * x) % m == 1)
            {
                return x;
            }
        }
        return 1;
    }
};

template <int Mod>
int niyuan(const ModInt<Mod> &x)
{
    return x.inverse().get_value();
}

int a[N];

void solve()
{
    int nn, kk;
    cin >> nn >> kk;
    ModInt<mod> n, k;
    n = nn, k = kk;
    ModInt<mod> sum = 0, cnt = 0;
    for (int i = 1; i <= nn; i++)
    {
        cin >> a[i];
        sum += a[i];
        if (a[i] != 0)
        {
            cnt += 1;
        }
    }
    ModInt<mod> fenmu = niyuan(n);
    ModInt<mod> p = cnt * fenmu;
    ModInt<mod> ans = sum * fenmu * k;
    for (int i = 1; i <= kk; i++)
    {
        ans = (ans + (k - i) * p * fenmu);
    }
    cout << ans << endl;
}

signed main()
{
    buff;
    int tt = 1;
    cin >> tt;
    while (tt--)
    {
        solve();
    }
    return 0;
}