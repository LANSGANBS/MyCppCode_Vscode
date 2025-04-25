# Capps Template Library

# 1, Math
## MathPackage
### Sieve
```cpp
class Sieve {
    std::vector<int> _mpf; // minimum prime factor
    std::vector<int> _primes;

    Sieve(std::size_t n) {
        init(n + 1);
    }

public:
    static Sieve &shared(std::size_t m = 128) {
        static Sieve instance(m);
        return instance;
    }

    void init(std::size_t n) {
        n = std::max(n, 2 * _mpf.size());
        if (_mpf.size() >= n) {
            return;
        }
        _mpf.assign(n, 0);
        _primes.clear();
        _primes.reserve(n / 10);
        for (auto i = 2U; i < n; i++) {
            if (_mpf[i] == 0) {
                _primes.push_back(i);
                _mpf[i] = i;
            }
            for (auto j = 0U; j < _primes.size() and i * _primes[j] < n; j++) {
                _mpf[i * _primes[j]] = _primes[j];
                if (_primes[j] == _mpf[i]) {
                    break;
                }
            }
        }
    }

    constexpr std::size_t size() const {
        return _mpf.size();
    }

    int mpf(std::size_t n) {
        if (n >= this->size()) {
            init(n);
        }
        return _mpf[n];
    }

    const std::vector<int> &primes() const {
        return _primes;
    }

    template <class T>
    std::vector<std::pair<T, int>> primeFactorize(T x) const {
        assert(1 <= x);
        const T n = static_cast<T>(this->size());
        std::vector<std::pair<T, int>> ps;
        auto process = [&](int d) {
            assert(d >= 2);
            int cnt = 0;
            while (x % d == 0) {
                x /= d;
                cnt += 1;
            }
            ps.emplace_back(d, cnt);
        };

        for (int d : _primes) {
            if (x < n or 1LL * d * d > x) {
                break;
            }

            if (x % d == 0) {
                process(d);
            }
        }

        for (int d = n; 1LL * d * d <= x and x >= n; d++) {
            if (x % d == 0) {
                process(d);
            }
        }

        if (x < n) {
            while (x > 1) {
                int d = _mpf[x];
                process(d);
            }
        }

        if (x > 1) {
            ps.emplace_back(x, 1);
        }

        return ps;
    }

    template <class T>
    std::vector<T> allFactors(T x) const {
        return allFactors(primeFactorize(x));
    }

    template <class T>
    static std::vector<T> allFactors(const std::vector<std::pair<T, int>> &ps) {
        std::vector<T> ds = {1};
        for (auto pr : ps) {
            T d = pr.first;
            int cnt = pr.second;
            int l = 0, r = ds.size();
            while (cnt--) {
                for (int k = l; k < r; k++) {
                    ds.push_back(ds[k] * d);
                }
                l = r, r = ds.size();
            }
        }

        return ds;
    }
};

auto siv = Sieve::shared(12345);
```

### ModuloInteger
```cpp
template <class T>
constexpr T power(T a, long long b) {
    T res = 1;
    for (; b; b /= 2, a *= a) {
        if (b & 1) {
            res *= a;
        }
    }
    return res;
}

template <class T, T P>
class ModuloInteger {
    T x;
    static T Mod;
    
    static constexpr int mult(int a, int b, const int Mod) {
        return 1LL * a * b % Mod;
    }

    static constexpr long long mult(long long a, long long b, const long long p) {
        long long res = a * b - static_cast<long long>(1.L * a * b / p) * p;
        res %= p;
        if (res < 0) {
            res += p;
        }
        return res;
    }

    constexpr T norm(T x) const {
		return (x < 0 ? x + getMod() : (x >= getMod() ? x - getMod() : x));
    }

public:
    typedef T ValueType;

    constexpr ModuloInteger() : x{} {}
    constexpr ModuloInteger(long long x) : x{norm(x % getMod())} {}

    static constexpr T getMod() {
		return (P > 0 ? P : Mod);
    }

    static constexpr void setMod(T Mod_) {
        Mod = Mod_;
    }

    constexpr T val() const {
        return x;
    }

    explicit constexpr operator T() const {
        return x;
    }

    constexpr ModuloInteger operator-() const {
        return ModuloInteger(getMod() - x);
    }

    constexpr ModuloInteger power(long long m) const {
        if (m < 0) {
            return ::power(inv(), -m);
        }
        return ::power((*this), m);
    }

    constexpr ModuloInteger inv() const {
        assert(x != 0);
        return this->power(getMod() - 2);
    }

    ModuloInteger &operator*=(ModuloInteger rhs) & {
        x = mult(x, rhs.x, getMod());
        return *this;
    }
    ModuloInteger &operator+=(ModuloInteger rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    ModuloInteger &operator-=(ModuloInteger rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    ModuloInteger &operator/=(ModuloInteger rhs) & {
        return *this *= rhs.inv();
    }

    friend constexpr ModuloInteger operator+(ModuloInteger lhs, ModuloInteger rhs) {
        return lhs += rhs;
    }

    friend constexpr ModuloInteger operator-(ModuloInteger lhs, ModuloInteger rhs) {
        return lhs -= rhs;
    }

    friend constexpr ModuloInteger operator*(ModuloInteger lhs, ModuloInteger rhs) {
        return lhs *= rhs;
    }

    friend constexpr ModuloInteger operator/(ModuloInteger lhs, ModuloInteger rhs) {
        return lhs /= rhs;
    }

    friend constexpr std::istream &operator>>(std::istream &is, ModuloInteger &a) {
        long long v;
        is >> v;
        a = ModuloInteger(v);
        return is;
    }

    friend constexpr std::ostream &operator<<(std::ostream &os, const ModuloInteger &a) {
        return os << a.val();
    }

    friend constexpr bool operator==(ModuloInteger lhs, ModuloInteger rhs) {
        return lhs.val() == rhs.val();
    }

    friend constexpr bool operator!=(ModuloInteger lhs, ModuloInteger rhs) {
        return lhs.val() != rhs.val();
    }
};

template <>
int ModuloInteger<int, 0>::Mod = 998244353;

template <>
long long ModuloInteger<long long, 0>::Mod = 4179340454199820289;

constexpr int P = 998244353;

using Z = ModuloInteger<std::remove_cv<decltype(P)>::type, P>;
```

```plain
V <= 1e9 : 1004535809
V <= 1e15 : 1337006139375617
V <= 4e18 : 4179340454199820289
```

### Combinatorics
```cpp
template <class T>
class Comb {
    int n;
    std::vector<T> _jc, _ijc;

    void _checkSize(int x) {
        if (x > n) {
            init(x * 2);
        }
    }

    constexpr Comb() 
        : n{0},
        _jc{1},
        _ijc{1} {}

    explicit Comb(int m) : Comb() {
        init(m);
    }
public:
    static Comb &shared(int m = 8) {
        static Comb instance(m);
        return instance;
    }

    void init(int m) {
        m = std::min(m, T::getMod() - 1);
        if (m <= n) {
            return;
        }

        _jc.resize(m + 1);
        _ijc.resize(m + 1);

        for (int i = n + 1; i <= m; i++) {
            _jc[i] = _jc[i - 1] * i;
        }
        _ijc.back() = _jc.back().inv();
        for (int i = m; i > n; i--) {
            _ijc[i - 1] = _ijc[i] * i;
        }

        n = m;
    }

    T jc(int x) {
        _checkSize(x);
        return _jc[x];
    }

    T ijc(int x) {
        _checkSize(x);
        return _ijc[x];
    }

    T A(int a, int b) {
        if (a < b or b < 0) {
            return 0;
        }
        _checkSize(a);
        return _jc[a] * _ijc[a - b];
    }
    T C(int a, int b) {
        if (a < b or b < 0) {
            return 0;
        }
        _checkSize(a);
        return _jc[a] * _ijc[a - b] * _ijc[b];
    }
};
auto comb = Comb<Z>::shared();
```

### FloatPointNumber
```cpp
template <class T>
class FloatPointNumber {
    static constexpr T EPS = 1E-12;
    static_assert(EPS >= 0, "EPS < 0");

    static constexpr int sgn(T x) {
        return x < -EPS ? -1 : x > EPS;
    }

    static int precision;
    static std::string inputStr;

    T x;
public:
    constexpr FloatPointNumber() : x{} {}
    constexpr FloatPointNumber(T x) : x{x} {}

    constexpr T val() const {
        return x;
    }

    constexpr int sgn() const {
        return sgn(x);
    }

    template <class G>
    constexpr G round() const {
        return G(x + 0.5 + EPS);
    }

    static constexpr void setprecision(int len) {
        precision = len;
    }

    // 四则运算

    FloatPointNumber &operator+=(FloatPointNumber a) & {
        x += a.x;
        return *this;
    }

    friend constexpr FloatPointNumber operator+(FloatPointNumber a, FloatPointNumber b) {
        return a += b;
    }

    constexpr FloatPointNumber operator-() const {
        return FloatPointNumber(-x);
    }

    FloatPointNumber &operator-=(FloatPointNumber a) & {
        x = x - a.x;
        return *this;
    }

    friend constexpr FloatPointNumber operator-(FloatPointNumber a, FloatPointNumber b) {
        return a -= b;
    }

    FloatPointNumber &operator*=(FloatPointNumber a) & {
        x *= a.x;
        return *this;
    }

    friend constexpr FloatPointNumber operator*(FloatPointNumber a, FloatPointNumber b) {
        return a *= b;
    }

    constexpr FloatPointNumber &operator/=(FloatPointNumber a) & {
        x /= (long double)a.x;
        return *this;
    }

    friend constexpr FloatPointNumber operator/(FloatPointNumber a, FloatPointNumber b) {
        return a /= b;
    }

    // 比较运算

    friend constexpr int operator<(FloatPointNumber a, FloatPointNumber b) {
        return sgn(a.x - b.x) < 0;
    }

    friend constexpr int operator<=(FloatPointNumber a, FloatPointNumber b) {
        return sgn(a.x - b.x) <= 0;
    }

    friend constexpr int operator>(FloatPointNumber a, FloatPointNumber b) {
        return sgn(a.x - b.x) > 0;
    }

    friend constexpr int operator>=(FloatPointNumber a, FloatPointNumber b) {
        return sgn(a.x - b.x) >= 0;
    }

    friend constexpr bool operator==(FloatPointNumber a, FloatPointNumber b) {
        return sgn(a.x - b.x) == 0;
    }

    friend constexpr bool operator!=(FloatPointNumber a, FloatPointNumber b) {
        return sgn(a.x - b.x) != 0;
    }

    // 输入输出

    friend constexpr std::istream &operator>>(std::istream &is, FloatPointNumber &a) {
        is >> inputStr;
        if (std::is_same<T, long double>::value) {
            a = FloatPointNumber(std::stold(inputStr));
        } else {
            a = FloatPointNumber(std::stod(inputStr));
        }
        return is;
    }

    friend constexpr std::ostream &operator<<(std::ostream &os, FloatPointNumber a) {
        return os << std::fixed << std::setprecision(precision) << a.val();
    }

    // 常数
    static constexpr FloatPointNumber PI = FloatPointNumber(acosl(-1));
};

template <class T>
std::string FloatPointNumber<T>::inputStr;

template <class T>
constexpr FloatPointNumber<T> std::abs(FloatPointNumber<T> x) {
    if (x.val() < 0) {
        x = -x;
    }
    return x;
}

template <class T>
constexpr FloatPointNumber<T> std::atan2(FloatPointNumber<T> x, FloatPointNumber<T> y) {
    return std::atan2(1.L * x.val(), 1.L * y.val());
}

template <class T>
constexpr FloatPointNumber<T> std::sin(FloatPointNumber<T> x) {
    return std::sin(1.L * x.val());
}

template <class T>
constexpr FloatPointNumber<T> std::cos(FloatPointNumber<T> x) {
    return std::cos(1.L * x.val());
}

template <class T>
constexpr FloatPointNumber<T> std::sqrt(FloatPointNumber<T> x) {
    return std::sqrt(1.L * x.val());
}

template <class T>
int FloatPointNumber<T>::precision = 6;

using Float = FloatPointNumber<double>;
```

### ExGcd 
```cpp
template <class T>
T exgcd(T a, T b, T &x, T &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    T g = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return g;
}
```



对于方程 $ ax+by=c $ , 调用 exgcd , 求出 $ x_0 $ 和 $ y_0 $ 使得 $ ax_0+by_0=\gcd(a, b) $



则在 $ \gcd(a,b)\mid c $ 的情况下有通解



$ x = x_0 \times \frac{c}{\gcd(a, b)}+k\times \frac{b}{\gcd(a, b)}
\\ \ 
\\  \\
y = y_0 \times \frac{c}{\gcd(a, b)}-k\times \frac{a}{\gcd(a, b)}
 $

## Polynomial
### FastFourierTransform
```cpp
template <class T, template <class G> class Complex>
class Polynomial : public std::vector<T> {
    using Comp = Complex<T>;

    static std::vector<Comp> w[2];
    static std::vector<int> r;

    static void init(int _log) {
        if (r.size() == (1U << _log)) {
            return;
        }

        int n = 1 << _log;
        r.assign(n, 0);
        for (int i = 1; i < n; i++) {
            r[i] = (r[i >> 1] >> 1) | ((i & 1) << (_log - 1));
        }

        w[0].assign(n, Comp());
        w[1].assign(n, Comp());

        const T PI = acosl(-1);
        for (int i = 0; i < n; i++) {
            auto th = PI * i / n;
            auto cth = std::cos(1.L * th);
            auto sth = std::sin(1.L * th);
            w[0][i] = Comp(cth, sth);
            w[1][i] = Comp(cth, -sth);
        }
    }

    static void fft(std::vector<Comp> &a, int op) {
        int n = a.size();
        init(std::__lg(n));
        for (int i = 0; i < n; i++) {
            if (i < r[i]) {
                std::swap(a[i], a[r[i]]);
            }
        }
        for (int mid = 1; mid < n; mid <<= 1) {
            const int d = n / mid;
            for (int R = mid << 1, j = 0; j < n; j += R) {
                for (int k = 0; k < mid; k++) {
                    Comp x = a[j + k];
                    Comp y = w[op][d * k] * a[j + mid + k];
                    a[j + k] = x + y;
                    a[j + mid + k] = x - y;
                }
            }
        }
    }

public:
    using std::vector<T>::vector;

    constexpr friend Polynomial operator*(const Polynomial &a, const Polynomial &b) {
        if (a.size() == 0 or b.size() == 0) {
            return Polynomial();
        }
        int n = a.size() + b.size() - 1;
        int _log = std::__lg(2 * n - 1);
        int s = 1 << _log;
        if (std::min(a.size(), b.size()) < 128) {
            Polynomial res(n);
            for (auto i = 0U; i < a.size(); i++) {
                for (auto j = 0U; j < b.size(); j++) {
                    res[i + j] += a[i] * b[j];
                }
            }
            return res;
        }

        std::vector<Comp> p(s), q(s);
        for (auto i = 0U; i < a.size(); i++) {
            p[i] = Comp(a[i], 0);
        }
        for (auto i = 0U; i < b.size(); i++) {
            q[i] = Comp(b[i], 0);
        }

        fft(p, 0);
        fft(q, 0);
        for (int i = 0; i < s; i++) {
            p[i] = p[i] * q[i];
        }
        fft(p, 1);

        Polynomial res(n);
        for (int i = 0; i < n; i++) {
            res[i] = p[i].real() / s; // 默认浮点数
        }
        return res;
    }
};

template <class T, template <class G> class Complex>
std::vector<Complex<T>> Polynomial<T, Complex>::w[2] = {};

template <class T, template <class G> class Complex>
std::vector<int> Polynomial<T, Complex>::r;

using Float = double;
using Poly = Polynomial<Float, std::complex>;
```

### NumberTheoreticTransform
```cpp
template <class T>
struct Polynomial : public std::vector<T> {

    static std::vector<T> w;
    static constexpr auto P = T::getMod();

    static void initW(int r) {
        if (static_cast<int>(w.size()) >= r) {
            return;
        }

        w.assign(r, 0);
        w[r >> 1] = 1;
        T s = ::power(T(3), (P - 1) / r);
        for (int i = r / 2 + 1; i < r; i++) {
            w[i] = w[i - 1] * s;
        }
        for (int i = r / 2 - 1; i > 0; i--) {
            w[i] = w[i * 2];
        }
    }

    constexpr friend void dft(Polynomial &a) {
        const int n = a.size();
        assert((n & (n - 1)) == 0);
        initW(n);

        for (int k = n >> 1; k; k >>= 1) {
            for (int i = 0; i < n; i += k << 1) {
                for (int j = 0; j < k; j++) {
                    T v = a[i + j + k];
                    a[i + j + k] = (a[i + j] - v) * w[k + j];
                    a[i + j] = a[i + j] + v;
                }
            }
        }
    }

    constexpr friend void idft(Polynomial &a) {
        const int n = a.size();
        assert((n & (n - 1)) == 0);
        initW(n);

        for (int k = 1; k < n; k <<= 1) {
            for (int i = 0; i < n; i += k << 1) {
                for (int j = 0; j < k; j++) {
                    T x = a[i + j];
                    T y = a[i + j + k] * w[j + k];
                    a[i + j + k] = x - y;
                    a[i + j] = x + y;
                }
            }
        }

        a *= P - (P - 1) / n;

        std::reverse(a.begin() + 1, a.end());
    }

public:
    using std::vector<T>::vector;

    constexpr Polynomial truncate(int k) const {
        Polynomial p = *this;
        p.resize(k);
        return p;
    }

    constexpr friend Polynomial operator+(const Polynomial &a, const Polynomial &b) {
        Polynomial p(std::max(a.size(), b.size()));
        for (auto i = 0U; i < a.size(); i++) {
            p[i] += a[i];
        }
        for (auto i = 0U; i < b.size(); i++) {
            p[i] += b[i];
        }
        return p;
    }

    constexpr friend Polynomial operator-(const Polynomial &a, const Polynomial &b) {
        Polynomial p(std::max(a.size(), b.size()));
        for (auto i = 0U; i < a.size(); i++) {
            p[i] += a[i];
        }
        for (auto i = 0U; i < b.size(); i++) {
            p[i] -= b[i];
        }
        return p;
    }

    constexpr friend Polynomial operator-(const Polynomial &a) {
        int n = a.size();
        Polynomial p(n);
        for (int i = 0; i < n; i++) {
            p[i] = -a[i];
        }
        return p;
    }

    constexpr friend Polynomial operator*(T a, Polynomial b) {
        for (auto i = 0U; i < b.size(); i++) {
            b[i] *= a;
        }
        return b;
    }
    constexpr friend Polynomial operator*(Polynomial a, T b) {
        for (auto i = 0U; i < a.size(); i++) {
            a[i] *= b;
        }
        return a;
    }

    constexpr friend Polynomial operator/(Polynomial a, T b) {
        b = b.inv();
        for (auto i = 0U; i < a.size(); i++) {
            a[i] *= b;
        }
        return a;
    }

    constexpr Polynomial mulxk(int k) const {
        assert(k >= 0);
        Polynomial b = (*this);
        b.insert(b.begin(), k, 0);
        return b;
    }

    constexpr Polynomial modxk(int k) const {
        assert(k > 0);
        return Polynomial(this->begin(), this->begin() + k);
    }
    
    constexpr Polynomial divxk(int k) const {
        assert(k >= 0);
        if (static_cast<int>(this->size()) <= k) {
            return Polynomial{};
        }

        return Polynomial(this->begin() + k, this->end());
    }

    constexpr T whenXis(T x) const {
        T ans = T{};
        for (int i = static_cast<int>(this->size()) - 1; i >= 0; i--) {
            ans = ans * x + this->at(i);
        }

        return ans;
    }

    Polynomial &operator+=(Polynomial b) {
        return (*this) = (*this) + b;
    }
    Polynomial &operator-=(Polynomial b) {
        return (*this) = (*this) - b;
    }
    Polynomial &operator*=(Polynomial b) {
        return (*this) = (*this) * b;
    }
    Polynomial &operator*=(T b) {
        return (*this) = (*this) * b;
    }
    Polynomial &operator/=(T b) {
        return (*this) = (*this) / b;
    }

    constexpr friend Polynomial operator*(const Polynomial &a, const Polynomial &b) {
        if (a.size() == 0 or b.size() == 0) {
            return Polynomial();
        }

        int n = a.size() + b.size() - 1;
        int s = 1 << std::__lg(2 * n - 1);

        if (((P - 1) & (s - 1)) != 0 or std::min(a.size(), b.size()) < 128) {
            Polynomial p(n);
            for (auto i = 0U; i < a.size(); i++) {
                for (auto j = 0U; j < b.size(); j++) {
                    p[i + j] += a[i] * b[j];
                }
            }

            return p;
        }

        Polynomial f = a.truncate(s);
        Polynomial g = b.truncate(s);

        dft(f), dft(g);

        for (int i = 0; i < s; i++) {
            f[i] *= g[i];
        }

        idft(f);

        return f.truncate(n);
    }

    constexpr Polynomial deriv() const {
        int n = this->size();
        if (n <= 1) {
            return Polynomial();
        }
        Polynomial p(n - 1);
        for (int i = 1; i < n; i++) {
            p[i - 1] = i * this->at(i);
        }
        return p;
    }

    constexpr Polynomial integr() const {
        int n = this->size();
        Polynomial p(n + 1);

        std::vector<T> _inv(n + 1);
        _inv[1] = 1;
        for (int i = 2; i <= n; i++) {
            _inv[i] = _inv[P % i] * (P - P / i);
        }

        for (int i = 0; i < n; ++i) {
            p[i + 1] = this->at(i) * _inv[i + 1];
        }
        return p;
    }

    // assert(this->at(0) != 0);
    constexpr Polynomial inv(int m = -1) const {
        const int n = this->size();
        m = m < 0 ? n : m;

        Polynomial p = Polynomial{this->at(0).inv()};
        p.reserve(4 * m);

        for (int k = 2; k / 2 < m; k <<= 1) {
            Polynomial q = Polynomial(this->begin(), this->begin() + std::min(k, n)).truncate(2 * k);
            p.resize(2 * k);

            dft(q), dft(p);
            for (int i = 0; i < 2 * k; i++) {
                p[i] = p[i] * (2 - p[i] * q[i]);
            }
            idft(p);

            p.resize(k);
        }

        return p.truncate(m);
    }

    constexpr Polynomial ln(int m = -1) const {
        m = m < 0 ? this->size() : m;

        return (deriv() * inv(m)).integr().truncate(m);
    }

    constexpr Polynomial exp(int m = -1) const {
        m = m < 0 ? this->size() : m;

        Polynomial p{1};

        int k = 1;
        while (k < m) {
            k <<= 1;
            p = (p * (Polynomial{1} - p.ln(k) + truncate(k))).truncate(k);
        }
        return p.truncate(m);
    }

    constexpr Polynomial power(long long k, int m, long long k2 = -1) const {
        if (0 < k and k <= (1 << 10)) {
            Polynomial p = (*this);
            Polynomial ans{1};
            for (; k; k >>= 1) {
                if (k & 1) {
                    ans *= p;
                    if (static_cast<int>(ans.size()) > m) {
                        ans.truncate(m);
                    }
                }

                p *= p;
                if (static_cast<int>(p.size()) > m) {
                    p.truncate(m);
                }
            }
            return ans.truncate(m);
        }

        k2 = k2 < 0 ? k : k2;
        unsigned int i = 0;
        while (i < this->size() and this->at(i) == T{}) {
            i++;
        }
        if (i == this->size() or k * i >= m) {
            return Polynomial(m, T{});
        }
        T v = this->at(i);
        Polynomial f = divxk(i) / v;
        return (f.ln(m - i * k) * k).exp(m - i * k).mulxk(i * k) * ::power(v, k2);
    }

    constexpr Polynomial sqrt(int m = -1) const {
        m = m < 0 ? this->size() : m;

        Polynomial p{1};

        int k = 1;
        constexpr T INV2 = T(1) / 2;
        while (k < m) {
            k <<= 1;
            p = (p + (truncate(k) * p.inv(k)).truncate(k)) * INV2;
        }

        return p.truncate(m);
    }

    friend constexpr std::istream &operator>>(std::istream &is, Polynomial &a) {
        int n = a.size();
        for (int i = 0; i < n; i++) {
            is >> a[i];
        }

        return is;
    }

    friend constexpr std::ostream &operator<<(std::ostream &os, const Polynomial &a) {
        int n = a.size();
        if (n >= 1) {
            os << a[0];
        }
        for (int i = 1; i < n; i++) {
            os << ' ' << a[i];
        }

        return os;
    }
};
template <class T>
std::vector<T> Polynomial<T>::w;

using Poly = Polynomial<Z>;
```

## Linear Algebra
### GaussianElimination
```cpp
template <class T>
std::string gauss(std::vector<std::vector<T>> &a) { // 传入增广矩阵
    const int n = a.size();
    const int m = a[0].size();
    assert(m > n);

    int c = 0, r = 0;
    for (; c < n; c++) { // c列r行，遍历列
        int tmp = r;
        for (int i = r; i < n; i++) { // 寻找列主元
            if (a[i][c] != 0) {
                tmp = i;
            }
        }
        if (a[tmp][c] == 0) { // 当前列全为0
            continue;
        }

        std::swap(a[tmp], a[r]); // 交换列主元

        for (int i = m - 1; i >= c; i--) { // 倒序处理
            a[r][i] /= a[r][c];
        }

        for (int i = r + 1; i < n; i++) {
            if (a[i][c] != 0) {
                for (int j = m - 1; j >= c; j--) {
                    a[i][j] -= a[r][j] * a[i][c];
                }
            }
        }

        r += 1;
    }

    if (r < n) {
        for (int i = r; i < n; i++) {
            if (a[i][n] != 0) {
                return "NoSolution";
            }
        }
        return "InfSolution";
    }

    for (int i = n - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            const T b = a[j][i];
            for (int k = n; k < m; k++) {
                a[j][k] -= b * a[i][k];
            }
            a[j][i] = 0;
        }
    }

    // 解会放在 a[i][n]  (0 <= i < n)
    return "OK";
}

template <class T>
struct MatrixUtil {
    std::string status;
    std::vector<std::vector<T>> inv;

    MatrixUtil(const std::vector<std::vector<T>> &a) {
        const int n = a.size();
        assert(static_cast<int>(a[0].size()) == n);

        std::vector<std::vector<T>> b(n, std::vector<T>(n * 2));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                b[i][j] = a[i][j];
            }
            b[i][i + n] = 1;
        }

        status = gauss(b);
        if (status == "OK") {
            inv = std::vector<std::vector<T>>(n, std::vector<T>(n));
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    inv[i][j] = b[i][j + n];
                }
            }
        }
    }
};
```

### LinearBasis
```cpp
template <class T>
struct LinearBasis {
    static constexpr int logn = std::__lg(std::numeric_limits<T>::max());
    std::array<T, logn + 1> b;
    int rank;
    bool canZero, dirt;

    LinearBasis() {
        clear();
    }

    void clear() {
        canZero = false;
        dirt = false;
        rank = 0;
        b.fill(0);
    }

    void insert(T x) {
        for (int i = logn; i >= 0; i--) {
            if (x >> i & 1) {
                if (b[i] == 0) {
                    b[i] = x;
                    rank += 1;
                    dirt = true;
                    return;
                }
                x ^= b[i];
            }
        }
        canZero = true;
    }

    // 询问线性基能不能异或出 x
    bool check(T x) {
        for (int i = logn; i >= 0; i--)
            if (x >> i & 1) {
                if (b[i] == 0) {
                    return false;
                }
                x ^= b[i];
            }
        return true;
    }

    T getMax() {
        T res = 0;
        for (int i = logn; i >= 0; i--) {
            res = std::max(res, res ^ b[i]);
        }
        return res;
    }

    T getMin() {
        if (canZero) {
            return 0;
        }
        int k = 0;
        while (k <= logn and b[k] == 0) {
            k += 1;
        }
        assert(k <= logn);
        return b[k];
    }

    void build() {
        if (dirt == false) {
            return;
        }
        for (int i = 1; i <= logn; i++) {
            for (int j = i - 1; j >= 0; j--) {
                if (b[i] >> j & 1) {
                    b[i] ^= b[j];
                }
            }
        }
        dirt = true;
    }

    T getDiffCount() {
        return (T(1) << rank) - !canZero;
    }

    // 查询有 K 个元素比它小的元素
    T find_by_order(T K) {
        if (canZero == false) {
            K += 1;
        }
        build();

        T res = 0;
        int j = 0;
        for (int i = 0; i <= logn; i++) {
            if (b[i] == 0) {
                continue;
            }

            if (K >> j & 1) {
                res ^= b[i];
            }
            j += 1;
        }

        return res;
    }
};
```

## Random Number Algorithm
### RandomNumber
```cpp
namespace __random {
    using u64 = unsigned long long;

    constexpr u64 chaos(u64 x) {
        return ((x ^ (x << 3)) ^ ((x ^ (x << 3)) >> 13)) ^
         (((x ^ (x << 3)) ^ ((x ^ (x << 3)) >> 13)) << 7);
    }

    constexpr u64 filter_string(u64 x, const char* str, size_t index) {
        return str[index] == '\0' ? x : filter_string(chaos(x ^ static_cast<u64>(str[index])), str, index + 1);
    }

    constexpr u64 generate_seed() {
        return filter_string(filter_string(filter_string(1128471 ^ __LINE__, __TIME__, 0), __TIMESTAMP__, 0), __FILE__, 0);
    };

    constexpr u64 seed = generate_seed();

    // __random float number
    template <class T>
    struct RandFloat {
        std::mt19937_64 myrand{seed};
        T operator()(T l, T r) {
            return std::uniform_real_distribution<T>(l, r)(myrand);
        }
    };
    using Float = double;
    __random::RandFloat<Float> randFloat;

    // __random integer number
    std::mt19937_64 rng(seed);
    // std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    // [l, r)
    template <class T>
    T randInt(T l, T r) {
        assert(l < r);
        return __random::rng() % (r - l) + l;
    }
};

using __random::randInt;
using __random::randFloat;
```

### MillerRabin
```cpp
/*
维基百科 :
n < 4e9, primes = [2, 7, 61]
n < 3e14, primes = [2, 3, 5, 7, 11, 13, 17]
n < 3e18, primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
n < 3e23, primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
*/
template <class T>
class MillerRabin {
    using R = ModuloInteger<T, 0>;
    const std::vector<int> primes;

    std::vector<int> getPrimes(int) {
        return std::vector<int>{2, 7, 61};
    }
    std::vector<int> getPrimes(long long) {
        return std::vector<int>{2, 3, 5, 7, 11, 13, 17, 19, 23};
    }

public:
    MillerRabin() : primes(getPrimes(T{})) {}

    constexpr bool operator()(T v) { // 判断v是不是质数
        if (v < 2 or (v != 2 and v % 2 == 0) or (v != 3 and v % 3 == 0)) {
            return false;
        }
        R::setMod(v);
        T s = v - 1;
        while ((s & 1) == 0) {
            s /= 2;
        }
        for (int p : primes) {
            if (v == p) {
                return true;
            }
            T t = s;
            R m = R(p).power(s);
            while (t + 1 != v and m - 1 != 0 and m + 1 != 0) {
                m = m * m;
                t *= 2;
            }

            if (m + 1 != 0 and (t & 1) == 0) {
                return false;
            }
        }
        return true;
    }
};
MillerRabin<long long> isPrime;
```

### PollardRho
```cpp
template <class T>
class PollardRho {
    using R = ModuloInteger<T, 0>;
    std::mt19937_64 rng{static_cast<unsigned long long>(std::chrono::steady_clock::now().time_since_epoch().count())};
    MillerRabin<T> rabin;

    static constexpr T gcd(T a, T b) {
        return b == 0 ? a : gcd(b, a % b);
    }

public:
    // 返回 n 的随机一个[2, n-1]内的因子, 或者判定是质数
    T findFactor(T n) {
        assert(n >= 2);
        if (n % 2 == 0) {
            return 2;
        }
        if (rabin(n)) {
            return n;
        }

        R::setMod(n);
        while (true) {
            T c = rng() % (n - 1) + 1;
            auto f = [&](R x) { return x * x + c; };
            R t = 0, r = 0;
            R p = 1, q;
            do {
                for (int i = 0; i < 128; i++) {
                    t = f(t);
                    r = f(f(r));
                    if (t == r or (q = p * std::abs(T(t) - T(r))) == 0) {
                        break;
                    }
                    p = q;
                }
                T d = gcd(T(p), n);
                if (d > 1) {
                    return d;
                }
            } while (t != r);
        }
    }

    std::vector<std::pair<T, int>> primeFactorize(T x) {
        assert(x >= 2);
        std::vector<T> ps;
        std::vector<T> stk = {x};
        while (!stk.empty()) {
            x = stk.back();
            stk.pop_back();

            T y = findFactor(x);
            if (x == y) {
                ps.push_back(x);
            } else {
                stk.push_back(y);
                stk.push_back(x / y);
            }
        }

        std::sort(begin(ps), end(ps));
        std::vector<std::pair<T, int>> pc;
        for (T p : ps) {
            if (pc.empty() or pc.back().first != p) {
                pc.emplace_back(p, 1);
            } else {
                pc.back().second += 1;
            }
        }

        return pc;
    }
};
PollardRho<long long> rho;
```

如果 $ n $ 是质数 $ (MillerRabbin判) $ 返回 $ n $

否则返回 $ n $ 的随机一个 $ [2,n-1] $ 的因子

复杂度理论 $ O(n^{\frac{1}{4}}\log n) $  但实际跑得快, 可以按 $ O(n^{\frac{1}{4}}) $ 算

# 2, Dynamic Programming
## MultipleBackpacks
```cpp
using VecGood = std::vector<std::array<int, 3>>;
std::vector<long long> multiBag(const VecGood &goods, int S) {
    // S 总背包大小
    std::vector<long long> f(S + 1, 0);
    std::vector<int> q(S + 2, 0);
    for (auto gd : goods) {
        int v = gd[0];
        int w = gd[1];
        int m = gd[2];
        // v价值, w体积, m数量
        auto calc = [&](int i, int j) {
            return f[j] + 1LL * (i - j) / w * v;
        };
        for (int up = S; up + w > S; up--) {
            int l = 1, r = 0, k = up;
            for (int x = up; x > 0; x -= w) {
                for (; k >= std::max(0LL, x - 1LL * m * w); k -= w) {
                    while (l <= r and calc(x, k) > calc(x, q[r])) {
                        r -= 1;
                    }
                    q[++r] = k;
                }
                f[x] = calc(x, q[l]);
                if (q[l] == x) {
                    l += 1;
                }
            }
        }
    }
    return f;
}
```

# 3, Sorting
## RadixSort
```cpp
template <int B, class T>
void radixSort(std::vector<T> &a) {
    if (a.empty()) {
        return;
    }
    const int mask = (1 << B) - 1, n = a.size();

    std::vector<T> b(n);
    std::vector<int> cnt(1 << B);

    T maxV = *std::max_element(begin(a), end(a));

    for (int i = 0; maxV; i += B, maxV >>= B) {
        std::fill(begin(cnt), end(cnt), 0);
        for (int j = 0; j < n; j++)
            cnt[a[j] >> i & mask] += 1;
        for (int j = 1; j < (1 << B); j++)
            cnt[j] += cnt[j - 1];
        for (int j = n - 1; j >= 0; j--)
            b[--cnt[a[j] >> i & mask]] = a[j];
        std::swap(a, b);
    }
}
```

## CountingSort
```cpp
// 按 a 的值返回下标排序
template <class T, class Key = int (*)(T)>
std::vector<int> countingSort(const std::vector<T> &a, Key key = [](T x) -> int { return x; }) {
    auto it = std::max_element(a.begin(), a.end(), [&](T p, T q) {
        return key(p) < key(q);
    });
    const int maxV = key(*it);
    std::vector<int> cnt(maxV + 1);

    for (auto i = 0U; i < a.size(); i++) {
        cnt[key(a[i])] += 1;
    }
    for (int i = 1; i <= maxV; i++) {
        cnt[i] += cnt[i - 1];
    }
    std::vector<int> res(a.size());
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; i--) {
        res[--cnt[key(a[i])]] = i;
    }

    return res;
}
```

## Discreter
```cpp
template <class T>
struct Discreter {
    std::vector<T> elementSet;

    // 待离散化的元素集合a
    Discreter(const std::vector<T> &a) : elementSet(a) {
        std::sort(begin(elementSet), end(elementSet));
        elementSet.erase(std::unique(begin(elementSet), end(elementSet)), end(elementSet));
    }

    std::vector<int> process(const std::vector<T> &a) const {
        std::vector<int> discRes(a.size());

        for (auto i = 0U; i < a.size(); i++) {
            discRes[i] = query(a[i]);
        }

        return discRes;
    }

    int query(const T &x) const {
        auto it = std::lower_bound(begin(elementSet), end(elementSet), x);
        // assert(it != end(elementSet) and *it == x);

        return it - begin(elementSet);
    }

    int queryUpperBound(const T &x) const {
        auto it = std::upper_bound(begin(elementSet), end(elementSet), x);

        return it - begin(elementSet);
    }

    T queryInv(int index) const {
        return elementSet[index];
    }

    int size() const {
        return elementSet.size();
    }
};
```

# 4, String
## String F4
### StringHash
```cpp
template <int D, const int *B, const int *P>
struct StringHash {
    std::vector<std::array<int, D>> h;

    template <class T>
    StringHash(const T &s) : h(s.size() + 1) {
        for (auto i = 0U; i < s.size(); i++) {
            for (int k = 0; k < D; k++) {
                h[i + 1][k] = (1LL * h[i][k] * B[k] + s[i] + 1) % P[k];
            }
        }
    }

    StringHash(const char *s) : StringHash(std::string(s)) {}

    // [l, r)
    std::array<int, D> get(int l, int r) {
        static std::vector<std::array<int, D>> spow(1);
        assert(l < r);

        if (static_cast<int>(spow.size()) < r - l + 1) {
            if (spow[0][0] == 0) {
                spow[0].fill(1);
            }
            int n = spow.size();
            spow.resize(r - l + 1);
            for (int i = n; i < static_cast<int>(spow.size()); i++) {
                for (int k = 0; k < D; k++) {
                    spow[i][k] = 1LL * spow[i - 1][k] * B[k] % P[k];
                }
            }
        }

        std::array<int, D> res = {};
        for (int k = 0; k < D; k++) {
            res[k] = h[r][k] - 1LL * h[l][k] * spow[r - l][k] % P[k];
            res[k] += (res[k] < 0 ? P[k] : 0);
        }
        return res;
    }
};

using u64 = unsigned long long;

namespace compileRandom {
    constexpr u64 chaos(u64 x) {
        return ((x ^ (x << 3)) ^ ((x ^ (x << 3)) >> 13)) ^
         (((x ^ (x << 3)) ^ ((x ^ (x << 3)) >> 13)) << 7);
    }

    constexpr u64 filter_string(u64 x, const char* str, size_t index) {
        return str[index] == '\0' ? x : filter_string(chaos(x ^ static_cast<u64>(str[index])), str, index + 1);
    }

    constexpr u64 generate_seed() {
        return filter_string(filter_string(filter_string(1128471 ^ __LINE__, __TIME__, 0), __TIMESTAMP__, 0), __FILE__, 0);
    };

    constexpr u64 seed = generate_seed();

    template <unsigned int T>
    struct Rng { static constexpr u64 value = chaos(Rng<T - 1>::value); };

    template <>
    struct Rng<0> { static constexpr u64 value = seed; };
}

constexpr int HashDimension = 2;

constexpr int __B[HashDimension] = {
    static_cast<int>(compileRandom::Rng<13>::value % 133 + 133),
    static_cast<int>(compileRandom::Rng<31>::value % 331 + 331)
};

constexpr int __P[HashDimension] = {
    static_cast<int>(1E9) + 21,
    static_cast<int>(1E9) + 33
};

using Hash = StringHash<HashDimension, __B, __P>;
```

### Kmp
```cpp
std::vector<int> kmp(const std::string &s) {
    int n = s.size();
    std::vector<int> link(n);
    for (int i = 1, j = 0; i < n; i++) {
        while (j and s[i] != s[j]) {
            j = link[j - 1];
        }   
        j += (s[i] == s[j]);
        link[i] = j;
    }
    return link;
}
```

### Z-Function
```cpp
std::vector<int> zFunction(const std::string &s) {
    int n = s.size();
    std::vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; i++) {
        if (i < r) {
            z[i] = std::min(z[i - l], r - i);
        }
        while (i + z[i] < n and s[i + z[i]] == s[z[i]]) {
            z[i]++;
        }
        if (i + z[i] > r) {
            l = i;
            r = i + z[i];
        }
    }
    return z;
}
```

### Manacher
```cpp
class Manacher {
    const int n;
    std::vector<int> r, f;
public:
    Manacher(const std::string &t)
        : n(t.size()), r(2 * t.size() + 3), f(2 * t.size() + 3) {
        std::string s = "[-";
        for (int i = 0; i < n; i++) {
            s += t[i];
            s += '-';
        }
        s.push_back(']');

        int mid = 1, far = 1;
        for (int i = 1; i < static_cast<int>(s.size()); i++) {
            r[i] = std::min(r[2 * mid - i], far - i);
            while (s[i + r[i]] == s[i - r[i]])
                r[i] += 1;
            if (far < i + r[i])
                mid = i, far = i + r[i];
            f[i + r[i] - 1] = std::max(f[i + r[i] - 1], r[i]);
        }
        for (int i = f.size() - 2; i; i--)
            f[i] = std::max(f[i], f[i + 1] - 1);
    }

    // 下标, 是否要以 +0.5 为中心
    int getPalinLenFromCenter(int i, int center) const {
        assert((!center and 0 <= i and i < n) or
               (center and 0 <= i and i < n - 1));

        return r[2 * (i + 1) + center] - 1;
    }

    int getPalinLenFromTail(int i) const {
        assert(0 <= i and i < n);
        return f[2 * (i + 1)];
    }
};
```

## Advanced Strings
### AcAutomaton
```cpp
template <int Z, char Base>
class AcAutomaton {
    void insert(int id, const std::string &s) {
        int p = 0;
        for (char c : s) {
            c -= Base;
            if (!son[p][c]) {
                son[p][c] = ++tot;
            }
            p = son[p][c];
        }
        // ID[p].push_back(id);
    }
    void build() {
        std::queue<int> q;
        for (int &y : son[0]) {
            if (y > 0) {
                q.push(y);
            }
        }
        while (!q.empty()) {
            int x = q.front();
            q.pop();

            int c = 0;
            for (int &y : son[x]) {
                if (y) {
                    link[y] = son[link[x]][c];
                    q.push(y);
                } else {
                    y = son[link[x]][c];
                }
                c += 1;
            }
        }
    }

public:
    std::vector<std::array<int, Z>> son;
    // std::vector<std::vector<int>> ID;
    std::vector<int> link;
    int tot = 0;

    AcAutomaton(const std::vector<std::string> &s) {
        int sizeSum = 0;
        for (auto t : s) {
            sizeSum += t.size();
        }
        son.resize(sizeSum + 1);
        // ID.resize(sizeSum + 1);
        link.resize(sizeSum + 1);

        for (auto i = 0U; i < s.size(); i++) {
            insert(i, s[i]);
        }
        build();
    }
};
```

### SuffixAutomaton
```cpp
template <int Z, char Base>
class Sam {
public:
    std::vector<std::array<int, Z>> son;
    std::vector<int> link, len;
    // std::vector<i64> cnt;
    int last;
    unsigned int tot;

    Sam(int n) : son(n * 2), link(n * 2), len(n * 2) {
        // cnt.assign(n * 2, 0);
        last = tot = 0;
        link[0] = -1;
    }

    Sam(const std::string &s) : Sam(s.size()) {
        for (auto i = 0U; i < s.size(); i++) {
            add(i, s[i]);
        }
    }

    void add(int id, char c) {
        c -= Base;
        int cur = ++tot;
        assert(tot < son.size());
        assign(id, cur, c);
        len[cur] = len[last] + 1;
        int v = last;
        while (v != -1 and !son[v][c]) {
            son[v][c] = cur;
            v = link[v];
        }

        if (v == -1) {
            link[cur] = 0;
        } else {
            int q = son[v][c];
            if (len[v] + 1 == len[q]) {
                link[cur] = q;
            } else {
                int clone = ++tot;
                assert(tot < son.size());
                len[clone] = len[v] + 1;
                son[clone] = son[q];
                link[clone] = link[q];

                while (v != -1 and son[v][c] == q) {
                    son[v][c] = clone;
                    v = link[v];
                }

                link[q] = link[cur] = clone;
            }
        }
        last = cur;
    }
private:
    void assign(int id, int cur, char c) {
    }
};
```

### ExSuffixAutomaton
```cpp
template <int Z, char Base>
class ExSam {
public:
    std::vector<std::array<int, Z>> son;
    std::vector<int> len, link;
    // std::vector<std::set<int>> ID;
    int tot = 0, last = 0;

    ExSam(const std::vector<std::string> &strs) {
        int sizeSum = 0;
        for (const auto &str : strs) {
            sizeSum += str.size();
        }

        len.resize(sizeSum * 2);
        link.resize(sizeSum * 2);
        link[0] = -1;
        son.resize(sizeSum * 2);
        // ID.resize(sizeSum * 2);

        for (auto i = 0U; i < strs.size(); i++) {
            last = 0;
            for (char c : strs[i]) {
                exadd(c, i);
            }
        }
    }
    void exadd(char c, int id) {
        c -= Base;
        if (son[last][c]) {
            int v = last, q = son[v][c];
            if (len[q] != len[v] + 1) {
                int cl = ++tot;
                len[cl] = len[v] + 1;
                link[cl] = link[q];
                son[cl] = son[q];

                while (v != -1 and son[v][c] == q)
                    son[v][c] = cl, v = link[v];
                link[q] = cl;
                q = cl;
            }
            int cur = last = q;
            assign(cur, id);
            return;
        }

        int cur = ++tot;
        len[cur] = len[last] + 1;
        assign(cur, id);

        int v = last;
        while (v != -1 and son[v][c] == 0)
            son[v][c] = cur, v = link[v];
        if (v == -1)
            link[cur] = 0;
        else {
            int q = son[v][c];
            if (len[q] == len[v] + 1)
                link[cur] = q;
            else {
                int cl = ++tot;
                len[cl] = len[v] + 1;
                link[cl] = link[q];
                son[cl] = son[q];

                while (v != -1 and son[v][c] == q)
                    son[v][c] = cl, v = link[v];
                link[q] = link[cur] = cl;
            }
        }
        last = cur;
    }

private:
    void assign(int cur, int id) {
        // ID[cur].insert(id);
    }
};
```

### PalindromeAutomanton
```cpp
template <int Z, char Base>
class Pam {
public:
    std::vector<std::array<int, Z>> son;
    std::vector<int> link, len, dep, cnt;
    std::string s;

    int cur = 0, tot = 1;

    Pam(int n)
        : son(n + 2), link(n + 2), len(n + 2),
          dep(n + 2), cnt(n + 2) {
        link[0] = 1;
        len[1] = -1;
    }

    Pam(const std::string &s) : Pam(s.size()) {
        for (auto i = 0U; i < s.size(); i++) {
            add(i, s[i]);
        }
    }

    int getLink(int x, int i) {
        while (i - len[x] - 1 < 0 or s[i - len[x] - 1] != s[i]) {
            x = link[x];
        }
        return x;
    }

    void add(int i, char c) {
        c -= Base;
        s.push_back(c);

        int v = getLink(cur, i);

        if (!son[v][c]) {
            link[++tot] = son[getLink(link[v], i)][c];
            son[v][c] = tot;
            len[tot] = len[v] + 2;
            dep[tot] = dep[link[tot]] + 1;
        }

        cur = son[v][c];
        assign(cur, i);
    }

    // Pam 的 linkTree 是 1 为根的
    std::vector<std::vector<int>> getLinkTree() const {
        std::vector<std::vector<int>> adj(tot + 1, std::vector<int>());
        for (int i = 0; i <= tot; i++) {
            if (i != 1) {
                adj[link[i]].push_back(i);
            }
        }
        return adj;
    }

private:
    void assign(int cur, int id) {
        cnt[cur] += 1;
    }
};
```

### SuffixArray
```cpp
class SuffixArray {
public:
    const int n;
    std::vector<int> sa, rk, h;

    template <class T>
    SuffixArray(const T &s)
        : n(s.size()), sa(n), rk(n), id(n), tmp(n) {

        std::iota(begin(id), end(id), 0);
        for (int i = 0; i < n; i++)
            rk[i] = s[i];

        countSort();

        for (int w = 1;; w <<= 1) {
            std::iota(begin(id), begin(id) + w + 1, n - w);
            for (int i = 0, p = w; i < n; i++)
                if (sa[i] >= w)
                    id[p++] = sa[i] - w;

            countSort();
            oldrk = rk;

            rk[sa[0]] = 0;
            for (int i = 1, p = 0; i < n; i++)
                rk[sa[i]] = equal(sa[i], sa[i - 1], w) ? p : ++p;

            if (rk[sa.back()] + 1 == n)
                break;
        }

        calcHeight(s);
    }

private:
    std::vector<int> oldrk, id, tmp, cnt;

    template <class T>
    void calcHeight(const T &s) {
        h.assign(n, 0);
        for (int i = 0, k = 0; i < n; i++) {
            if (rk[i] == 0)
                continue;
            k -= bool(k);
            while (s[i + k] == s[sa[rk[i] - 1] + k])
                k += 1;
            h[rk[i]] = k;
        }
    }

    // 计数排序
    void countSort() {
        int m = *std::max_element(begin(rk), end(rk));
        cnt.assign(m + 1, 0);
        for (int i = 0; i < n; i++)
            cnt[tmp[i] = rk[id[i]]] += 1;
        for (auto i = 1U; i < cnt.size(); i++)
            cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; i--)
            sa[--cnt[tmp[i]]] = id[i];
    }

    bool equal(int x, int y, int w) {
        int rkx = (x + w < n ? oldrk[x + w] : -1);
        int rky = (y + w < n ? oldrk[y + w] : -1);
        return oldrk[x] == oldrk[y] and rkx == rky;
    }

    /**
     * sa[i] -> 第i小的后缀
     * rk[i] -> 后缀i的排名
     */
};
```



# 5, DataStructures
## BaseDataStructures
### HashTable
```cpp
using u64 = unsigned long long;

template <class T, int Mod>
class HashTable {
    struct Iterator {
        u64 *pKey;
        T *pVal;

        Iterator(u64 *pKey, T *pVal) : pKey{pKey}, pVal{pVal} {}

        Iterator operator++() {
            ++pKey;
            ++pVal;

            return *this;
        }

        bool operator!=(Iterator it) const {
            return pKey != it.pKey;
        }

        std::pair<u64, T> operator*() {
            return std::make_pair(*pKey, *pVal);
        }
    };

    Iterator begin() const {
        return Iterator(to + 1, val + 1);
    }

    Iterator end() const {
        return Iterator(to + tot + 1, val + tot + 1);
    }

    int hd[Mod], nt[Mod * 2], tot = 0;
    u64 to[Mod * 2];
    T val[Mod * 2];
public:
    void clear() {
        for (int i = 1; i <= tot; i++) {
            hd[to[i] % Mod] = 0;
        }
        tot = 0;
    }

    T operator()(u64 x) const {
        int u = x % Mod;
        for (int i = hd[u]; i > 0; i = nt[i]) {
            if (to[i] == x) {
                return val[i];
            }
        }
        return T{};
    }

    T &operator[](u64 x) {
        int u = x % Mod;
        for (int i = hd[u]; i > 0; i = nt[i]) {
            if (to[i] == x) {
                return val[i];
            }
        }
        to[++tot] = x;
        nt[tot] = hd[u];
        hd[u] = tot;
        return val[tot] = T{};
    }
};

HashTable<int, int(1e6) + 3> mp;
```

### SparseTable
```cpp
template <class T, class Cmp = std::less<T>>
class RMQ {
    const int n;
    const int logn;

    const Cmp cmp = Cmp();
    std::vector<std::vector<T>> jump;
public:
    RMQ(const std::vector<T> &a)
        : n(a.size()), logn{std::__lg(n)}, jump(logn + 1) {

        jump[0] = a;

        for (int j = 1; j <= logn; j++) {
            jump[j].resize(n - (1 << j) + 1);
        }

        for (int j = 0; j < logn; j++) {
            const int limit = n - (1 << (j + 1));
            for (int i = 0; i <= limit; i++) {
                jump[j + 1][i] = std::min(jump[j][i], jump[j][i + (1 << j)], cmp);
            }
        }
    }

    // [l, r)
    constexpr T operator()(int l, int r) const {
        assert(l < r and r <= n);
        int log = std::__lg(r - l);
        return std::min(jump[log][l], jump[log][r - (1 << log)], cmp);
    }
};
```

### DeletableHeap
```cpp
template <class T, class Cmp = std::less<T>>
class DeletableHeap {
    std::priority_queue<T, std::vector<T>, Cmp> q_push, q_erase;
    // Heap = q_push - q_erase

    void check() {
        while (!q_erase.empty() and q_push.top() == q_erase.top()) {
            q_push.pop(), q_erase.pop();
        }
    }

public:
    void push(T x) {
        q_push.push(x);
    }

    void erase(T x) {
        q_erase.push(x);
    }

    T top() {
        check();
        return q_push.top();
    }

    void pop() {
        check();
        q_push.pop();
    }

    int size() {
        int q_push_size = q_push.size();
        int q_erase_size = q_erase.size();
        assert(q_push_size >= q_erase_size);
        return q_push_size - q_erase_size;
    }
};
```

### DisjointSetUnion
```cpp
class DSU {
    std::vector<int> f;
public:
    std::vector<int> size;

    DSU(int n) {
        init(n);
    }

    void init(int n) {
        f.assign(n, 0);
        std::iota(f.begin(), f.end(), 0);

        size.assign(n, 1);
    }

    int find(int x) {
        while (x != f[x])
            x = f[x] = f[f[x]];
        return x;
    }

    void Union(int x, int y) {
        x = find(x), y = find(y);
        if (x == y)
            return;

        // if (size[x] < size[y]) {
        //     std::swap(x, y);
        // }

        size[x] += size[y];
        f[y] = x;
    }
};
```

## TreeDataStructures
### LazySegmentTree
```cpp
template <class T, class Tag>
class LazySegT {
    int n;
    std::vector<T> info;
    std::vector<Tag> tag;

    void up(int i) {
        info[i] = info[i * 2] + info[i * 2 + 1];
    }
    void apply(int i, const Tag &v) {
        info[i].apply(v);
        tag[i].apply(v);
    }
    void down(int i) {
        apply(i * 2, tag[i]);
        apply(i * 2 + 1, tag[i]);
        tag[i] = Tag();
    }

    void modify(int i, int l, int r, int x, const T &v) {
        if (r - l == 1) {
            info[i] = v;
            return;
        }
        int mid = (l + r) / 2;
        down(i);
        if (x < mid) {
            modify(i * 2, l, mid, x, v);
        } else {
            modify(i * 2 + 1, mid, r, x, v);
        }
        up(i);
    }

    T rangeQuery(int i, int l, int r, int tl, int tr) {

        if (tl <= l and r <= tr)
            return info[i];

        down(i);
        int mid = (l + r) / 2;

        return (tl < mid ? rangeQuery(i * 2, l, mid, tl, tr) : T()) +
               (mid < tr ? rangeQuery(i * 2 + 1, mid, r, tl, tr) : T());
    }

    void rangeModify(int i, int l, int r, int tl, int tr, const Tag &v) {

        if (tl <= l and r <= tr) {
            apply(i, v);
            return;
        }
        down(i);
        int mid = (l + r) / 2;

        if (tl < mid)
            rangeModify(i * 2, l, mid, tl, tr, v);
        if (mid < tr)
            rangeModify(i * 2 + 1, mid, r, tl, tr, v);
        up(i);
    }

    template <class F>
    int findFirst(int i, int l, int r, int tl, int tr, F pred) {
        if (l >= tr || r <= tl || !pred(info[i])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int mid = (l + r) / 2;
        down(i);
        int res = findFirst(i * 2, l, mid, tl, tr, pred);
        if (res == -1) {
            res = findFirst(i * 2 + 1, mid, r, tl, tr, pred);
        }
        return res;
    }

    template <class F>
    int findLast(int i, int l, int r, int tl, int tr, F pred) {
        if (l >= tr || r <= tl || !pred(info[i])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int mid = (l + r) / 2;
        down(i);
        int res = findLast(i * 2 + 1, mid, r, tl, tr, pred);
        if (res == -1) {
            res = findLast(i * 2, l, mid, tl, tr, pred);
        }
        return res;
    }
public:
    LazySegT(int n, T v = T())
        : LazySegT(std::vector<T>(n, v)) {}

    template <class G>
    LazySegT(const std::vector<G> &a) : n(a.size()) {
        info.assign(4 << std::__lg(n), T());
        tag.assign(4 << std::__lg(n), Tag());
        std::function<void(int, int, int)> build = [&](int i, int l, int r) {
            if (r - l == 1) {
                info[i] = T(a[l]);
                return;
            }
            int mid = (l + r) / 2;
            build(i * 2, l, mid);
            build(i * 2 + 1, mid, r);
            up(i);
        };
        build(1, 0, n);
    }

    // 单点赋值
    void modify(int i, const T &v) {
        modify(1, 0, n, i, v);
    }

    // 区间查询 [l, r)
    T rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }

    // 区间修改 [l, r)
    void rangeModify(int l, int r, const Tag &v) {
        return rangeModify(1, 0, n, l, r, v);
    }

    // 区间左边第一个满足条件的下标
    template <class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }

    // 区间右边第一个满足条件的下标
    template <class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 0, n, l, r, pred);
    }
};

struct Tag {
    int add;

    explicit Tag(int add = 0)
        : add(add){}

    void apply(const Tag &tag) {
        if (tag.add) {
            add += tag.add;
        }
    }
};

struct Node {
    i64 val;
    int len;

    explicit Node(i64 val = 0, int len = 1)
        : val(val), len(len) {}

    void apply(const Tag &tag) {
        if (tag.add) {
            val += 1LL * tag.add * len;
        }
    }

    Node operator+(const Node &a) {
        return Node(val + a.val, len + a.len);
    }
};
```

### PresidentTree
```cpp
template <class T>
class PresidentTree {
    using NodeIndex = int;
    struct Node {
        int val;
        NodeIndex l, r;

        Node(int val = 0)
            : val{val}, l{0}, r{0} {}
    };

    std::vector<Node> t; // memory pool

    std::vector<NodeIndex> root;
    const T Start, Last;

    NodeIndex newNode(int val = 0) {
        t.push_back(val);
        return (int)t.size() - 1;
    }

    void up(NodeIndex i) {
        t[i].val = t[t[i].l].val + t[t[i].r].val;
    }

    void modify(NodeIndex &p, T l, T r, T x) {
        if (p == 0) {
            p = newNode();
        }
        if (r - l == 1) {
            t[p].val++;
            return;
        }

        T mid = (0LL + l + r) / 2;

        if (x < mid)
            modify(t[p].l, l, mid, x);
        else
            modify(t[p].r, mid, r, x);
        up(p);
    }
    NodeIndex merge(NodeIndex x, NodeIndex y, T l, T r) {
        if (!x or !y)
            return (x ? x : y);

        // 每次把 x 修改
        if (r - l == 1) {
            t[x].val += t[y].val;
            return x;
        }

        T mid = (0LL + l + r) / 2;
        t[x].l = merge(t[x].l, t[y].l, l, mid);
        t[x].r = merge(t[x].r, t[y].r, mid, r);
        return up(x), x;
    }

    constexpr int getRange(NodeIndex x, NodeIndex y, T l, T r, T tl, T tr) {
        if (tl <= l and r <= tr) {
            return t[y].val - t[x].val;
        }
        T mid = (0LL + l + r) / 2;
        return (tl < mid ? getRange(t[x].l, t[y].l, l, mid, tl, tr) : 0) +
               (mid < tr ? getRange(t[x].r, t[y].r, mid, r, tl, tr) : 0);
    }

    constexpr T getKth(NodeIndex x, NodeIndex y, T l, T r, int k) {
        if (r - l == 1)
            return l;
        T mid = (0LL + l + r) / 2;
        int L = t[t[y].l].val - t[t[x].l].val;

        return (L >= k ? getKth(t[x].l, t[y].l, l, mid, k)
                       : getKth(t[x].r, t[y].r, mid, r, k - L));
    }

public:
    PresidentTree(const std::vector<T> &a, T min, T max)
        : t(1), root(a.size() + 1), Start(min), Last(max + 1) {

        t.reserve(a.size() * std::__lg(a.size() * 2));

        root[0] = newNode();
        for (auto i = 1U; i <= a.size(); i++) {
            if (t.capacity() <= t.size() + 64) {
                t.reserve(std::max(2 * t.capacity(), t.capacity() + 64));
            }
            modify(root[i], Start, Last, a[i - 1]);
            root[i] = merge(root[i], root[i - 1], Start, Last);
        }
    }
    // [l, r), [tl, tr)
    constexpr int getRange(int l, int r, T tl, T tr) {
        return getRange(root[l], root[r], Start, Last, tl, tr);
    }
    // [l, r)
    constexpr T getKth(int l, int r, int k) {
        return getKth(root[l], root[r], Start, Last, k);
    }
};
```

### FenwickTree
```cpp
template <class T, class Cmp = std::greater<T>>
struct Max {
    const Cmp cmp = Cmp();
    constexpr T operator()(const T &a, const T &b) const {
        return std::min(a, b, cmp);
    }
};

template <class T, class Merge = std::plus<T>>
class Fenwick {
    const int n;
    std::vector<T> t;
    const Merge merge;
public:
    Fenwick(int n) : n{n}, t(n + 1), merge(Merge()) {}

    // O(n) build Fenwick
    Fenwick(const std::vector<T> &a) : Fenwick(a.size()) {
        for (int i = 1; i <= n; i++) {
            t[i] = merge(t[i], a[i - 1]);
            int j = i + (i & -i);
            if (j <= n) {
                t[j] = merge(t[j], t[i]);
            }
        }
    }

    void modify(int i, const T &x) {
        for (i += 1; i <= n; i += i & -i) {
            t[i] = merge(t[i], x);
        }
    }

    T posQuery(int i) const {
        T res = T{};
        for (i += 1; i > 0; i -= i & -i) {
            res = merge(res, t[i]);
        }
        return res;
    }

    // [l, r)
    // T rangeQuery(int l, int r) {
    //     return posQuery(r - 1) - posQuery(l - 1);
    // }

    // 合并起来 <= k 的最长前缀
    int select(const T &k) const {
        int x = 0;
        T cur{};
        for (int i = 1 << std::__lg(n); i > 0; i /= 2) {
            if (x + i <= n and merge(cur, t[x + i]) <= k) {
                x += i;
                cur = merge(cur, t[x]);
            }
        }
        return x - 1;
    }
};
```

# 6, Tree Theory
## TreePre
```cpp
template <class T>
class TreePre {
    static constexpr int endPoint(int x) {
        return x;
    }
    template <class G>
    static constexpr int endPoint(const std::pair<int, G> &pr) {
        return pr.first;
    }

    void dfs1(int x, int f) {
        fa[x] = f;

        for (auto p : adj[x]) {
            int y = endPoint(p);
            if (y != f) {
                dep[y] = dep[x] + 1;
                dfs1(y, x);

                size[x] += size[y];
                if (big[x] == -1 or size[y] > size[big[x]]) {
                    big[x] = y;
                }
            }
        }
    }
    void dfs2(int x, int top) {
        dfn[x] = cur++;
        idfn[dfn[x]] = x;
        tp[x] = top;
        if (big[x] != -1) {
            dfs2(big[x], top);
        }

        for (auto p : adj[x]) {
            int y = endPoint(p);
            if (y != big[x] and y != fa[x]) {
                dfs2(y, y);
            }
        }
    }
    const std::vector<std::vector<T>> &adj;

    const int n;
    const int root;

public:
    std::vector<int> big, size, tp, dep, fa, dfn, idfn;
    // dfn begin from 0
    int cur = 0;

    TreePre(const std::vector<std::vector<T>> &g, int root)
        : adj(g), n(g.size()), root{root}, big(n, -1), size(n, 1),
          tp(n), dep(n), fa(n), dfn(n), idfn(n) {
        // dep begin from 0
        // dep[root] = 0;
        dfs1(root, -1);
        dfs2(root, root);
    }

    int getLca(int x, int y) const {
        while (tp[x] != tp[y]) {
            if (dep[tp[x]] > dep[tp[y]]) {
                x = fa[tp[x]];
            } else {
                y = fa[tp[y]];
            }
        }
        return (dep[x] < dep[y] ? x : y);
    }

    int dist(int x, int y) const {
        int lca = getLca(x, y);
        return dep[x] + dep[y] - 2 * dep[lca];
    }

    // x→y路径剖分的dfn号区间[l, r], l > r 说明这是上升段
    std::vector<std::pair<int, int>> getRoad(int x, int y) const {
        int lca = getLca(x, y);
        std::vector<std::pair<int, int>> vec1, vec2;
        while (tp[x] != tp[lca]) {
            vec1.push_back({dfn[x], dfn[tp[x]]});
            x = fa[tp[x]];
        }

        if (x != lca) {
            vec1.push_back({dfn[x], dfn[lca] + 1});
        }

        vec1.push_back({dfn[lca], dfn[lca]});

        while (tp[y] != tp[lca]) {
            vec2.push_back({dfn[tp[y]], dfn[y]});
            y = fa[tp[y]];
        }

        if (y != lca) {
            vec2.push_back({dfn[lca] + 1, dfn[y]});
            y = fa[tp[y]];
        }

        vec1.insert(vec1.end(), vec2.rbegin(), vec2.rend());
        return vec1;
    }

    int kthAncester(int x, int k) const {
        if (dep[x] - dep[root] < k) {
            return -1;
        }

        int d = dep[x] - k;

        while (dep[tp[x]] > d) {
            x = fa[tp[x]];
        }

        return idfn[dfn[x] - dep[x] + d];
    }

    // x is y's ancester
    bool isAncester(int x, int y) const {
        return dfn[x] <= dfn[y] and dfn[y] < dfn[x] + size[x];
    }
};
```

## VirtualTreePre
```cpp
// T:原树边, G:虚树边
template <class T, class G>
class VirtualTreePre {
    const std::function<bool(int, int)> dfnCmp = [&](int x, int y) {
        return pre.dfn[x] < pre.dfn[y];
    };

public:
    const TreePre<T> pre;
    std::vector<std::vector<G>> vt;

    VirtualTreePre(const std::vector<std::vector<T>> &adj, int root)
        : pre(adj, root), vt(adj.size()) {}

    // 虚树存在vt, 返回vt根节点
    int build(std::vector<int> a) {
        std::sort(a.begin(), a.end(), dfnCmp);
        for (int i = (int)a.size() - 1; i; i--) {
            a.push_back(pre.getLca(a[i - 1], a[i]));
        }
        std::sort(a.begin(), a.end(), dfnCmp);
        a.erase(std::unique(a.begin(), a.end()), a.end());

        for (int x : a) {
            vt[x].clear();
        }

        for (auto i = 1U; i < a.size(); i++) {
            int lca = pre.getLca(a[i - 1], a[i]);
            vt[lca].emplace_back(a[i], pre.dep[a[i]] - pre.dep[lca]);
        }

        return a.front();
    }
};
```

## CentroidDecomposition
```cpp
class CentroidDecomposition {
    const std::vector<std::vector<int>> &adj;
    std::vector<bool> vis;
    std::vector<int> size;

    void dfsSize(int x, int fa) {
        size[x] = 1;
        for (int y : adj[x]) {
            if (y != fa and !vis[y]) {
                dfsSize(y, x);
                size[x] += size[y];
            }
        }
    }

    int getRoot(int x, int fa, int m) {
        for (int y : adj[x]) {
            if (y != fa and !vis[y] and 2 * size[y] > m) {
                return getRoot(y, x, m);
            }
        }
        return x;
    }

    void build(int x) {
        vis[x] = true;
        dfsOrder.push_back(x);

        for (int y : adj[x]) {
            if (!vis[y]) {
                dfsSize(y, -1);
                y = getRoot(y, -1, size[y]);

                // cdt[x].push_back(y);
                build(y);
            }
        }
    }

public:
    // std::vector<std::vector<int>> cdt;
    std::vector<int> dfsOrder;
    int root;
    CentroidDecomposition(const std::vector<std::vector<int>> &g)
        : adj(g),
          vis(g.size()),
          //   cdt(g.size()),
          size(g.size()) {
        dfsSize(0, -1);
        root = getRoot(0, -1, size[0]);
        build(root);
    }
};
```

# 7, Graph Theory
## GridPre
```cpp
// 网格单位元素, 网格一行Vector, 网格输出邻接表单位元素
template <class Element, class Vector, class AdjEle>
class GridPre {
    const int h, w, n;

    static constexpr int D = 4;
    static constexpr int dx[D] = {0, 0, -1, 1};
    static constexpr int dy[D] = {-1, 1, 0, 0};

public:
    std::vector<std::vector<AdjEle>> adj;
    // std::vector<Element> a;
    std::vector<int> reachableVertexSet;

    GridPre(const std::vector<Vector> &grid, const std::set<Element> &reachableSet)
        : h{grid.size()}, w{grid[0].size()}, n{h * w},
          //   a(n),
          adj(n) {

        reachableVertexSet.reserve(n);

        auto ok = [&](int x, int y) {
            return 0 <= x and x < h and 0 <= y and y < w and reachableSet.count(grid[x][y]);
        };

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                if (reachableSet.count(grid[i][j]) == 0) {
                    continue;
                }

                int x = i * w + j;
                // a[x] = grid[i][j];

                reachableVertexSet.push_back(x);

                for (int k = 0; k < D; k++) {
                    int next_i = i + dx[k];
                    int next_j = j + dy[k];

                    if (ok(next_i, next_j)) {
                        adj[x].emplace_back(next_i * w + next_j);
                    }
                }
            }
        }
    }
};
```

## Dijkstra
```cpp
template <class T, class G>
class Dijkstra {
    const std::vector<std::vector<std::pair<int, T>>> &adj;
    std::vector<std::vector<G>> dis;

    std::vector<G> get(int s) {
        std::vector<G> dis(adj.size(), std::numeric_limits<G>::max() / 2);

        using Pair = std::pair<G, int>;
        std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> q;

        dis[s] = G();
        q.push({dis[s], s});

        while (!q.empty()) {
            auto _pr = q.top();
            G nearDist = _pr.first;
            int x = _pr.second;
            q.pop();

            if (nearDist > dis[x]){
                continue;
            }

            for (auto pr : adj[x]) {
                int y = pr.first;
                T w = pr.second;
                if (dis[y] > dis[x] + w) {
                    dis[y] = dis[x] + w;
                    q.push({dis[y], y});
                }
            }
        }
        return dis;
    }

public:
    Dijkstra(const std::vector<std::vector<std::pair<int, T>>> &g)
        : adj(g), dis(g.size()) {}

    G operator()(int x, int y) {
        if (dis[x].empty())
            dis[x] = get(x);
        return dis[x][y];
    }
};
```

## RingTree
```cpp
class RingTree {
public:
    std::vector<std::vector<int>> g;
    std::vector<int> ring;
    // DSU dsu;

    RingTree(const std::vector<std::vector<int>> &adj)
        : g(adj.size())
    //   dsu(adj.size())
    {
        const int n = adj.size();

        std::vector<int> deg(n, 0);
        std::queue<int> q;
        for (int x = 0; x < n; x++) {
            deg[x] = adj[x].size();
            if (deg[x] == 1) {
                q.push(x);
            }
        }

        while (!q.empty()) {
            int x = q.front();
            q.pop();

            for (int y : adj[x]) {
                deg[y] -= 1;
                if (deg[y] == 1) {
                    q.push(y);
                }
                if (deg[y] >= 1) {
                    g[y].push_back(x);
                    // dsu.Union(y, x);
                }
            }
        }

        int x = 0;
        while (deg[x] < 2) {
            x += 1;
        }
        assert(x < n and deg[x] == 2);

        do {
            ring.push_back(x);
            int z = -1;
            for (int y : adj[x]) {
                if (deg[y] == 2 and y != x) {
                    deg[x] -= 1;
                    deg[y] -= 1;
                    z = y;
                    break;
                }
            }
            x = z;
        } while (x >= 0);
    }
};
```

## TopSort
```cpp
class TopSort {
    static constexpr int endPoint(int x) {
        return x;
    }
    template <class G>
    static constexpr int endPoint(const std::pair<int, G> &pr) {
        return pr.first;
    }

public:
    template <class T>
    std::vector<int> operator()(const std::vector<T> &adj) const {
        int n = adj.size();
        std::vector<int> ind(n);
        for (int x = 0; x < n; x++) {
            for (auto p : adj[x]) {
                ind[endPoint(p)] += 1;
            }
        }

        std::vector<int> q;
        for (int x = 0; x < n; x++) {
            if (ind[x] == 0) {
                q.push_back(x);
            }
        }

        std::vector<int> res;
        while (!q.empty()) {
            int x = q.back();
            res.push_back(x);
            q.pop_back();

            for (auto p : adj[x]) {
                int y = endPoint(p);
                ind[y] -= 1;
                if (ind[y] == 0) {
                    q.push_back(y);
                }
            }
        }

        return res;
    }
};
const TopSort topSort;
```

## Connectivity
### StronglyConnectedComponent
```cpp
class SCC {
    const std::vector<std::vector<int>> &adj;
    std::vector<int> q; // stack
    int r = 0, cur = 0;

    void dfs(int x) {
        dfn[x] = low[x] = cur++;
        q[++r] = x;

        for (int y : adj[x]) {
            if (dfn[y] == -1) {
                dfs(y);
                low[x] = std::min(low[x], low[y]);
            } else if (bel[y] == -1) {
                low[x] = std::min(low[x], dfn[y]);
            }
        }

        if (dfn[x] == low[x]) {
            int y;
            do {
                y = q[r--];
                bel[y] = cntBlock;
            } while (y != x);
            cntBlock += 1;
        }
    }

public:
    // original graph
    std::vector<int> dfn, low, bel;

    // shrinking graph
    std::vector<std::vector<int>> g;
    int cntBlock = 0;

    SCC(const std::vector<std::vector<int>> &adj)
        : adj(adj), dfn(adj.size(), -1), low(adj.size()), bel(adj.size(), -1) {
        int n = adj.size();
        q.assign(n + 1, 0);

        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                dfs(i);
            }
        }

        g.resize(cntBlock);
        for (int x = 0; x < n; x++) {
            for (int y : adj[x]) {
                if (bel[x] == bel[y])
                    continue;
                g[bel[x]].push_back(bel[y]);
            }
        }

        // for (int x = 0; x < cntBlock; x++) {
        //     std::sort(begin(g[x]), end(g[x]));
        //     g[x].erase(std::unique(begin(g[x]), end(g[x])), end(g[x]));
        // }
    }
};
```

### TwoSat
```cpp
class TwoSat {
    const int n;
    std::vector<std::vector<int>> adj;
public:
    std::vector<bool> ans;
    TwoSat(int n) : n(n), adj(2 * n), ans(n) {}
    // x * 2         为 x 设 0
    // x * 2 + 1     为 x 设 1

    // 如果 x 取 f 那么 y 要取 g
    void add(int x, bool f, int y, bool g) {
        x *= 2, y *= 2;
        adj[x + f].push_back(y + g);
        adj[y + !g].push_back(x + !f);
    }

    void assign(int x, bool f) {
        x *= 2;
        adj[x + !f].push_back(x + f);
    }

    bool work() {
        SCC scc(adj);
        const auto &bel = scc.bel;

        for (int i = 0; i < n; i++) {
            if (bel[2 * i] == bel[2 * i + 1]) {
                return false;
            }   
            ans[i] = bel[2 * i] > bel[2 * i + 1];
        }

        return true;
    }
};
```

### VertexBiconnectedComponent
```cpp
class VertexBC {
    const int n;
    const std::vector<std::vector<int>> &adj;
    std::stack<int, std::vector<int>> q;
    int cur = 0, sqid;

    void dfs(int x, int root) {
        dfn[x] = low[x] = cur++;
        q.push(x);

        for (int y : adj[x]) {
            if (dfn[y] == -1) {
                dfs(y, root);
                low[x] = std::min(low[x], low[y]);

                if (low[y] == dfn[x]) {
                    csqt.push_back({});
                    for (int z = -1; z != y; q.pop()) {
                        z = q.top();
                        csqt[z].push_back(sqid);
                        csqt[sqid].push_back(z);
                    }
                    csqt[x].push_back(sqid);
                    csqt[sqid].push_back(x);
                    sqid += 1;
                }
            } else {
                low[x] = std::min(low[x], dfn[y]);
            }
        }
    }

public:
    // original graph
    std::vector<int> dfn, low;
    std::vector<std::vector<int>> csqt; // 圆方树
    int componentNum = 0;

    VertexBC(const std::vector<std::vector<int>> &adj)
        : n(adj.size()),
          adj(adj),
          sqid{n},
          dfn(n, -1),
          low(n),
          csqt(n)
           {

        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                componentNum += 1;
                dfs(i, i);
                q.pop();
            }
        }
    }
};
```

### EdgeBiconnectedComponent
```cpp
class EdgeBC {
    const std::vector<std::vector<int>> &adj;
    std::vector<int> q; // stack
    int r = 0, cur = 0;

    void dfs(int x, int fa) {
        dfn[x] = low[x] = cur++;
        q[++r] = x;

        for (int y : adj[x]) {
            if (y == fa) {
                fa = ~fa;
                continue;
            }
            if (dfn[y] == -1) {
                dfs(y, x);
                low[x] = std::min(low[x], low[y]);
            } else {
                low[x] = std::min(low[x], dfn[y]);
            }
        }

        if (dfn[x] == low[x]) {
            int y;
            do {
                y = q[r--];
                bel[y] = cntBlock;
            } while (y != x);
            cntBlock += 1;
        }
    }

public:
    // original graph
    std::vector<int> dfn, low, bel, cutDeg;

    // shrinking graph
    std::vector<std::vector<int>> g;
    int cntBlock = 0, componentNum = 0;

    EdgeBC(const std::vector<std::vector<int>> &adj)
        : adj(adj), dfn(adj.size(), -1), low(adj.size()), bel(adj.size(), -1), cutDeg(adj.size()) {
        int n = adj.size();
        q.assign(n + 1, 0);

        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                componentNum += 1;
                dfs(i, -1);
            }
        }

        g.resize(cntBlock);
        for (int x = 0; x < n; x++) {
            for (int y : adj[x]) {
                if (bel[x] == bel[y])
                    continue;
                g[bel[x]].push_back(bel[y]);
            }
        }
    }
};
```

## Flow
### MaxFlow
```cpp
template <class T>
class Flow {
    std::vector<int> cur, dep;
    bool bfs(int s, int t) {
        dep.assign(n, -1);
        std::queue<int> q;
        dep[s] = 0;

        q.push(s);
        while (!q.empty()) {
            const int u = q.front();
            q.pop();

            for (int i : g[u]) {
                auto pr = adj[i];
                int v = pr.first;
                T c = pr.second;

                if (c > 0 and dep[v] == -1) {
                    dep[v] = dep[u] + 1;
                    if (v == t)
                        return true;
                    q.push(v);
                }
            }
        }

        return false;
    }

    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        T res = f;
        for (int &i = cur[u]; i < static_cast<int>(g[u].size()); i++) {
            const int j = g[u][i];
            auto pr = adj[j];
            int v = pr.first;
            T c = pr.second;

            if (c > 0 and dep[v] == dep[u] + 1) {
                T out = dfs(v, t, std::min(res, c));
                adj[j].second -= out;
                adj[j ^ 1].second += out;

                res -= out;
                if (res == 0) {
                    return f;
                }
            }
        }
        return f - res;
    }

public:
    int n;
    std::vector<std::pair<int, T>> adj;
    std::vector<std::vector<int>> g;
    static constexpr T Inf = std::numeric_limits<T>::max();

    Flow(int m) : n(m), g(m) {}

    int newNode() {
        n += 1;
        g.push_back({});
        return n - 1;
    }

    int add(int u, int v, T c, T c_ = 0) {
        int eid = adj.size();
        g[u].push_back(adj.size());
        adj.emplace_back(v, c);
        g[v].push_back(adj.size());
        adj.emplace_back(u, c_);

        return eid;
    }

    T work(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, Inf);
        }
        return ans;
    }

    std::vector<std::pair<int, int>> getCuts(int s, int t) {
        std::string vis(n, 't');
        vis[s] = 's';

        std::queue<int> q;
        q.push(s);

        while (!q.empty()) {
            int x = q.front();
            q.pop();
            for (int i : g[x]) {
                auto pr = adj[i];
                int y = pr.first;
                T w = pr.second;
                if (vis[y] == 's' or w == 0) {
                    continue;
                }
                vis[y] = 's';
                q.push(y);
            }
        }

        std::vector<std::pair<int, int>> cuts;
        for (int x = 0; x < n; x++) {
            if (vis[x] == 't') {
                continue;
            }

            for (int i : g[x]) {
                auto pr = adj[i];
                int y = pr.first;
                T w = pr.second;
                if (vis[y] == 't') {
                    assert(w == 0);
                    cuts.emplace_back(x, y);
                }
            }
        }

        return cuts;
    }
};
```

### MinCostFlow
```cpp
template <class T, class F>
class MCFGraph {
    struct Edge {
        int y;
        T c;
        F f;
        Edge(int y, T c, F f) : y{y}, c{c}, f{f} {}
    };

    bool dijkstra(int s, int t) {
        dis.assign(n, std::numeric_limits<F>::max());
        pre.assign(n, -1);
        using Pair = std::pair<F, int>;
        std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> q;
        dis[s] = 0;
        q.emplace(0, s);

        while (!q.empty()) {
            auto _pr = q.top();
            F nearDist = _pr.first;
            int x = _pr.second;
            q.pop();

            if (dis[x] < nearDist)
                continue;
            for (int i : g[x]) {
                int y = adj[i].y;
                T c = adj[i].c;
                F f = adj[i].f;
                if (c > 0 and dis[y] > nearDist + h[x] - h[y] + f) {
                    dis[y] = nearDist + h[x] - h[y] + f;
                    pre[y] = i;
                    q.emplace(dis[y], y);
                }
            }
        }
        return dis[t] != std::numeric_limits<F>::max();
    }
public:
    const int n;
    std::vector<Edge> adj;
    std::vector<std::vector<int>> g;
    std::vector<F> h, dis;
    std::vector<int> pre;
    
    MCFGraph(int n) : n(n), g(n) {}
    void add(int x, int y, T c, F f) {
        g[x].push_back(adj.size());
        adj.emplace_back(y, c, f);
        g[y].push_back(adj.size());
        adj.emplace_back(x, 0, -f);
    }
    std::pair<T, F> work(int s, int t) {
        T flow = 0;
        F cost = 0;
        h.assign(n, 0);
        while (dijkstra(s, t)) {
            for (int i = 0; i < n; ++i)
                h[i] += dis[i];
            T aug = std::numeric_limits<T>::max();
            for (int i = t; i != s; i = adj[pre[i] ^ 1].y)
                aug = std::min(aug, adj[pre[i]].c);
            for (int i = t; i != s; i = adj[pre[i] ^ 1].y) {
                adj[pre[i]].c -= aug;
                adj[pre[i] ^ 1].c += aug;
            }
            flow += aug;
            cost += F(aug) * h[t];
        }
        return std::pair<T, F>(flow, cost);
    }
};
```

# 8, Calculate Geometry
## Point
```cpp
template <class T, class G>
struct BaseVector2 {
    T x, y;
    constexpr BaseVector2() : BaseVector2(T{}, T{}) {}
    constexpr BaseVector2(T x, T y) : x{x}, y{y} {}

    // 运算
    constexpr BaseVector2 operator+(BaseVector2 a) const {
        return BaseVector2(x + a.x, y + a.y);
    }
    constexpr BaseVector2 operator-(BaseVector2 a) const {
        return BaseVector2(x - a.x, y - a.y);
    }
    constexpr BaseVector2 operator-() const {
        return BaseVector2(-x, -y);
    }
    constexpr G operator*(BaseVector2 a) const {
        return G(x) * a.x + G(y) * a.y;
    }
    constexpr G operator%(BaseVector2 a) const {
        return G(x) * a.y - G(y) * a.x;
    }
    constexpr BaseVector2 rotate() const {
        return BaseVector2(-y, x);
    }
    template <class F>
    constexpr BaseVector2 rotate(F theta) const {
        BaseVector2 b(std::cos(theta), std::sin(theta));
        return BaseVector2(x * b.x - y * b.y,
                           x * b.y + y * b.x);
    }
    constexpr friend BaseVector2 operator*(const T &a, BaseVector2 b) {
        return BaseVector2(a * b.x, a * b.y);
    }

    // 比较
    constexpr bool operator<(BaseVector2 a) const {
        if (x == a.x) {
            return y < a.y;
        }
        return x < a.x;
    }

    constexpr bool operator>(BaseVector2 a) const {
        if (x == a.x) {
            return y > a.y;
        }
        return x > a.x;
    }

    constexpr bool operator<=(BaseVector2 a) const {
        if (x == a.x) {
            return y <= a.y;
        }
        return x <= a.x;
    }

    constexpr bool operator>=(BaseVector2 a) const {
        if (x == a.x) {
            return y >= a.y;
        }
        return x >= a.x;
    }

    constexpr bool operator==(BaseVector2 a) const {
        return x == a.x and y == a.y;
    }

    constexpr bool operator!=(BaseVector2 a) const {
        return x != a.x or y != a.y;
    }

    // 输入输出
    friend std::istream &operator>>(std::istream &in, BaseVector2 &p) {
        return in >> p.x >> p.y;
    }
    friend std::ostream &operator<<(std::ostream &ot, BaseVector2 p) {
        return ot << '(' << p.x << ", " << p.y << ')';
    }
};

template <class T, class G>
G dis2(BaseVector2<T, G> a, BaseVector2<T, G> b = BaseVector2<T, G>(0, 0)) {
    BaseVector2<T, G> p = a - b;
    return p * p;
}
template <class T, class G>
auto dis(BaseVector2<T, G> a, BaseVector2<T, G> b = BaseVector2<T, G>(0, 0)) -> decltype(std::sqrt(G(1))) {
    return std::sqrt(dis2(a, b));
}

template <class T, class G>
int sgn(BaseVector2<T, G> p) {
    if (p.x < 0 or p.x == 0 and p.y > 0) {
        return 1;
    } else
        return 0;
    // 以41象限为先
}

template <class Vector>
bool polarCmp(Vector p, Vector q) {
    if (sgn(p) == sgn(q)) {
        if (p % q == 0) {
            return dis2(p) < dis2(q);
        } else {
            return p % q > 0;
        }
    } else {
        return sgn(p) < sgn(q);
    }
}

using Point = BaseVector2<int, long long>;
using Vector = Point;
using PS = std::vector<Point>;
```

## Line
```cpp
// a on Seg b-c
bool onSeg(Point a, Point b, Point c) {
    b = b - a;
    c = c - a;
    return b % c == 0 and b * c <= 0;
}

// a disTo Line b-c
Float disToLine(Point a, Point b, Point c) {
    Point v1 = b - c, v2 = a - c;
    return std::abs(v1 % v2) / std::sqrt(v1 * v1);
}

// a disTo Seg b-c
Float disToSeg(Point a, Point b, Point c) {
    if ((a - b) * (c - b) <= 0 or (a - c) * (b - c) <= 0) {
        return std::min(dis(a, b), dis(a, c));
    }
    return disToLine(a, b, c);
}

// a project to Line b-c
Point foot(Point a, Point b, Point c) {
    Point u = a - b, v = c - b;
    return b + (u * v) / (v * v) * v;
}

// a symmetry to Line b-c
Point symmetry(Point a, Point b, Point c) {
    Point ft = foot(a, b, c);
    return 2 * ft - a;
}

// Line a-b cross Line c-d
Point cross(Point a, Point b, Point c, Point d) {
    Point v = c - d;
    Float sa = v % (a - d), sb = (b - d) % v;
    return sa / (sa + sb) * (b - a) + a;
}

// a-b 和 a-c 的夹角 signed
Float getAngle(Point a, Point b, Point c) {
    auto v1 = b - a, v2 = c - a;
    return std::atan2(v1 % v2, v1 * v2); // ∠bac
}

// 对四个不同的点判断四点共圆
// d在abc外接圆外return 1, 内return -1
int inCircle(Point a, Point b, Point c, Point d) {
    const Float PI = acosl(-1);
    Float ag1 = getAngle(a, b, c), ag2 = getAngle(d, c, b);
    auto sgn = [](Float x) { return x < 0 ? -1 : (x > 0); };
    if (sgn(ag1) == sgn(ag2)) {
        return sgn(PI - std::abs(ag1 + ag2));
    } else {
        return sgn(std::abs(ag1) - std::abs(ag2));
    }
}
```

## Polygon
###  About ConvexHull 
```cpp
PS convex(PS ps) {
    std::sort(ps.begin(), ps.end());
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());

    const int n = ps.size();
    if (n <= 1) {
        return ps;
    }

    PS hull(n + 1);
    int k = -1;

    for (int i = 0; i < n; i++) {
        while (k >= 1 and (hull[k] - hull[k - 1]) % (ps[i] - hull[k]) <= 0) {
            k -= 1;
        }
        hull[++k] = ps[i];
    }

    for (int i = n - 2, m = k + 1; i >= 0; i--) {
        while (k >= m and (hull[k] - hull[k - 1]) % (ps[i] - hull[k]) <= 0) {
            k -= 1;
        }
        hull[++k] = ps[i];
    }

    assert(k >= 2);
    return hull.resize(k), hull;
}

template<class T>
T area(const PS &ps) {
    T res = 0;
    for (auto i = 0U; i < ps.size(); i++) {
        int j = (i + 1) % ps.size();
        res += ps[i] % ps[j];
    }
    return res / 2;
}

PS minkowski(const PS &a, const PS &b) {
    int n = a.size(), m = b.size();
    if (!n or !m)
        return PS();

    PS ps(1, a[0] + b[0]);
    int ap = 0, bp = 0;
    while (ap < n and bp < m) {
        auto va = a[(ap + 1) % n] - a[ap];
        auto vb = b[(bp + 1) % m] - b[bp];
        auto res = va % vb;
        if (res > 0)
            ps.push_back(ps.back() + va), ap++;
        if (res < 0)
            ps.push_back(ps.back() + vb), bp++;
        if (res == 0)
            ps.push_back(ps.back() + va + vb), ap++, bp++;
    }
    while (ap < n) {
        auto va = a[(ap + 1) % n] - a[ap];
        ps.push_back(ps.back() + va), ap++;
    }
    while (bp < m) {
        auto vb = b[(bp + 1) % m] - b[bp];
        ps.push_back(ps.back() + vb), bp++;
    }
    return ps.pop_back(), ps;
}
```

### Point-Line-Polygon
```cpp
// #include "onSeg"

class InPoly {
    // Seg c-d is Cross Seg a-b
    bool isCross(Point a, Point b, Point c, Point d) const {
        Vector ab = b - a, cd = d - c;
        if (ab % cd == 0)
            return 0; // 共线,寄
        int r1 = sgn(ab % (c - a)), r2 = sgn(ab % (d - a));
        int g1 = sgn(cd % (a - c)), g2 = sgn(cd % (b - c));
        return !(r1 * r2 > 0 or r1 + r2 == -1 or g1 * g2 > 0);
        // 真相交或者 c-d 贴着 a-b 左边
    }

public:
    // a IN/OUT/ON Polygon
    std::string operator()(Point a, const PS &ps) const {
        int res = 0;
        Vector b = {1 << 30, a.y};
        for (int i = 0; i < ps.size(); i++) {
            int j = (i + 1) % ps.size();
            if (onSeg(a, ps[i], ps[j]))
                return "ON";
            res += isCross(a, b, ps[i], ps[j]);
        }
        return (res % 2 ? "IN" : "OUT");
    }
};
const InPoly inPoly;

// a IN/OUT/ON Convex
std::string inConvex(Point a, const PS &ps) {
    if (a == ps[0])
        return "ON";
    if (ps.size() <= 1)
        return "OUT";
    if (ps.size() == 2)
        return onSeg(a, ps[0], ps[1]) ? "ON" : "OUT";
    auto v = a - ps[0];
    if ((ps[1] - ps[0]) % v < 0 or (ps.back() - ps[0]) % v > 0)
        return "OUT";
    int l = 1, r = ps.size() - 1;
    while (l + 1 < r) {
        auto mid = l + r >> 1;
        auto res = (ps[mid] - ps[0]) % v;
        if (res == 0)
            return (ps[mid] == a ? "ON" : (onSeg(a, ps[0], ps[mid]) ? "IN" : "OUT"));
        (res > 0 ? l : r) = mid;
    }
    auto res = (ps[r] - ps[l]) % (a - ps[l]);
    if (res == 0 or onSeg(a, ps[0], ps[l]) or onSeg(a, ps[0], ps[r]))
        return "ON";
    return (res > 0 ? "IN" : "OUT");
}

PS cutPoly(const PS &ps, Point a, Point b) {
    // 返回多边形 ps 在有向直线 a->b 左边的部分
    PS v;
    auto c = b - a;
    for (int i = 0; i < ps.size(); i++) {
        int j = (i + 1) % ps.size();
        auto cr1 = c % (ps[i] - a), cr2 = c % (ps[j] - a);
        if (cr1 >= 0)
            v.push_back(ps[i]);
        if (sgn(cr1) * sgn(cr2) == -1)
            v.push_back(cross(a, b, ps[i], ps[j]));
    }
    return v;
}


// find point of tangency
class BinarySearchOnConvex {
    PS vs;
public:
    constexpr BinarySearchOnConvex(const PS &ps) : vs(ps.size()) {
        int n = size(ps);
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            vs[i] = ps[j] - ps[i];
        }
    }

    int operator()(Point a) const {
        int res = std::lower_bound(begin(vs), end(vs), a, polarCmp) - begin(vs);
        return res % size(vs);
    }
};
```



# 9, TemplateFetcher
```python
# 配置表，将可配置的参数集中放在这里，方便后续修改
CONFIG = {
    # [[可配置]] 语雀文档的链接
    "yuque_url": "https://www.yuque.com/capps/ctl/heb78p9y3xvyrpz9",
    # [[可配置]] 目标目录，用于存放生成的 .code-snippets 文件
    "target_dir": "/mnt/c/Users/9600x/AppData/Roaming/Code/User/snippets/",
    # 自定义 code-snippets 配置数组
    "custom_snippets": [
        {
            "prefix": "_T_",
            "content": """\
#include <bits/stdc++.h>

using i64 = long long;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    $0

    return 0;
}
"""
        },
        {
            "prefix": "_T_case",
            "content": """\
#include <bits/stdc++.h>

using i64 = long long;

void solve() {
    $0

    return;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int t;
    std::cin >> t;

    while (t--) {
        solve();
    }

    return 0;
}
"""
        },
    ],
    
    # 中间产物目录名
    "_download_dir": "download",
    "_acm_templates_dir": "download/ACMtemplates",
    "_templates_dir": "templates",
    "_snippets_dir": "snippets",
}

import os
import json
import re
import subprocess
import shutil


# 检查并安装必要的工具
def check_and_install_tools():
    try:
        # 检查 npm 是否安装
        subprocess.run(['npm', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print("npm is installed.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("npm is not installed. Starting to install node/npm...")
        try:
            # 清楚原有依赖
            subprocess.run(['sudo', 'apt', 'remove', '--purge', 'nodejs', 'npm'], check=True)
            subprocess.run(['sudo', 'apt', 'autoremove'], check=True)
            # 更新软件包列表
            subprocess.run(['sudo', 'apt', 'update'], check=True)
            # 清理软件包缓存
            subprocess.run(['sudo', 'apt', 'clean'], check=True)
            # 修复依赖问题
            subprocess.run(['sudo', 'apt', 'install', '-f'], check=True)
            # 安装 Node.js 和 npm
            subprocess.run(['sudo', 'apt', 'install', 'build-essential', 'nodejs', '-y'], check=True)
            print("npm has been successfully installed.")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred during the installation: {e}")
            raise
    try:
        # 检查 yuque - dl 是否安装
        # https://github.com/gxr404/yuque-dl
        subprocess.run(['yuque-dl', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print("yuque-dl is installed.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("yuque-dl is not installed. Installing globally using npm...")
        try:
            # npm i -g yuque-dl // 安装yuque-dl
            subprocess.run(['npm', 'i', '-g', 'yuque-dl'], check=True)
            print("yuque-dl installed successfully.")
        except subprocess.CalledProcessError:
            print("Error occurred while installing yuque-dl.")
            raise


# 从语雀下载文档
def download_yuque_doc():
    try:
        subprocess.run(['yuque-dl', CONFIG["yuque_url"]], check=True)
    except subprocess.CalledProcessError:
        print("Error occurred while downloading Yuque document.")
        raise


# 从 Markdown 文件提取代码到 templates 目录
def extract_code_to_templates():
    download_dir = CONFIG["_download_dir"]
    acm_templates_dir = CONFIG["_acm_templates_dir"]
    templates_dir = CONFIG["_templates_dir"]

    # 检查并创建 templates 目录
    if not os.path.exists(templates_dir):
        os.makedirs(templates_dir)

    # 检查 ACMtemplates 目录是否存在
    if os.path.exists(acm_templates_dir):
        md_file_path = os.path.join(acm_templates_dir, "CappsTemplateLibrary.md")
        if os.path.exists(md_file_path):
            try:
                # 读取 Markdown 文件
                with open(md_file_path, 'r', encoding='utf-8') as file:
                    md_content = file.read()
                # 定义正则表达式模式，用于匹配 ```cpp 到 ``` 之间的代码块
                pattern = r'(\S+)\s+```cpp(.*?)```'
                matches = re.findall(pattern, md_content, re.DOTALL)
                # 遍历匹配结果
                for match in matches:
                    # 获取原始文件名
                    original_filename = match[0]
                    # 过滤非法字符，只保留字母、数字和下划线
                    valid_chars = ''.join(e for e in original_filename if e.isalnum() or e == '_')
                    # 拼接有效的文件名
                    filename = os.path.join(templates_dir, valid_chars + '.cpp')
                    # 获取代码内容
                    code = match[1].strip()
                    # 将代码写入对应的 .cpp 文件
                    with open(filename, 'w', encoding='utf-8') as cpp_file:
                        cpp_file.write(code)
                    # print(f"Code has been written to {filename}")
            except Exception as e:
                print(f"Error occurred while extracting code: {e}")
                raise
        else:
            print("The CappsTemplateLibrary.md file does not exist.")
            raise FileNotFoundError("The CappsTemplateLibrary.md file does not exist.")
    else:
        print("The ACMtemplates directory does not exist.")
        raise FileNotFoundError("The ACMtemplates directory does not exist.")


# 生成 code-snippets 文件
def generate_snippets():
    templates_dir = CONFIG["_templates_dir"]
    snippets_dir = CONFIG["_snippets_dir"]

    # 检查 snippets 目录是否存在，如果不存在则创建
    if not os.path.exists(snippets_dir):
        os.makedirs(snippets_dir)

    # 用于存储已使用的前缀
    used_prefixes = set()

    # 处理从 yuque 拉下来的文件生成 snippets
    for filename in os.listdir(templates_dir):
        if filename.endswith('.cpp'):
            # 获取文件名（不包含扩展名）
            base_name = os.path.splitext(filename)[0]
            prefix = f"_T_{base_name}"
            if prefix in used_prefixes:
                print(f"[[!error!]]Prefix {prefix} is already in use. Skipping...")
                continue
            used_prefixes.add(prefix)
            # 拼接 .cpp 文件的完整路径
            cpp_file_path = os.path.join(templates_dir, filename)
            # 拼接 .code - snippets 文件的完整路径
            snippets_file_path = os.path.join(snippets_dir, f'{base_name}.code-snippets')

            try:
                # 读取 .cpp 文件的内容，指定字符编码为 utf - 8
                with open(cpp_file_path, 'r', encoding='utf-8') as cpp_file:
                    cpp_content = cpp_file.readlines()

                # 处理 .cpp 文件内容，对 " 和 \ 进行转义
                escaped_content = []
                for line in cpp_content:
                    # 对双引号和反斜杠进行转义
                    line = line.replace('\\', '\\\\').replace('"', '\\"')
                    escaped_content.append(line.rstrip())

                # 构建 JSON 数据
                snippet = {
                    f"Print to console": {
                        "scope": "cpp",
                        "prefix": prefix,
                        "body": escaped_content + ["$0"],
                        "description": "Log output to console"
                    }
                }

                # 将 JSON 数据写入 .code - snippets 文件，指定字符编码为 utf - 8
                with open(snippets_file_path, 'w', encoding='utf-8') as snippets_file:
                    json.dump(snippet, snippets_file, indent=4, ensure_ascii=False)
                # print(f"Successfully generated {snippets_file_path}")

            except UnicodeDecodeError:
                print(f"Encoding error occurred while reading {cpp_file_path}. "
                      "Please ensure the file is in UTF-8 encoding.")
                raise
            except Exception as e:
                print(f"Error occurred while processing {cpp_file_path}: {e}")
                raise

    # 处理自定义的 code-snippets
    for custom_snippet in CONFIG["custom_snippets"]:
        prefix = custom_snippet["prefix"]
        if prefix in used_prefixes:
            print(f"[[!error!]]Prefix {prefix} is already in use. Skipping...")
            continue
        used_prefixes.add(prefix)
        content = custom_snippet["content"]
        escaped_content = []
        for line in content.splitlines():
            line = line.replace('\\', '\\\\').replace('"', '\\"')
            escaped_content.append(line.rstrip())

        # 构建 JSON 数据
        snippet = {
            "Print to console": {
                "scope": "cpp",
                "prefix": prefix,
                "body": escaped_content,
                "description": "Log output to console"
            }
        }

        # 拼接 .code - snippets 文件的完整路径
        snippets_file_path = os.path.join(snippets_dir, f'{prefix.replace("_T_", "")}.code-snippets')
        try:
            # 将 JSON 数据写入 .code - snippets 文件，指定字符编码为 utf - 8
            with open(snippets_file_path, 'w', encoding='utf-8') as snippets_file:
                json.dump(snippet, snippets_file, indent=4, ensure_ascii=False)
            print(f"Successfully generated {snippets_file_path}")
        except Exception as e:
            print(f"Error occurred while creating {snippets_file_path}: {e}")
            raise


# 将 snippets 目录下的文件移动到目标目录
def move_snippets_to_target():
    snippets_dir = CONFIG["_snippets_dir"]
    target_dir = CONFIG["target_dir"]
    # 检查源目录是否存在
    if not os.path.exists(snippets_dir):
        print("The source directory snippets does not exist.")
        raise FileNotFoundError("The source directory snippets does not exist.")

    # 检查目标目录是否存在，如果不存在则创建
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # 移动所有.code - snippets 文件到目标目录，允许覆盖
    for file in os.listdir(snippets_dir):
        if file.endswith('.code-snippets'):
            source_file = os.path.join(snippets_dir, file)
            target_file = os.path.join(target_dir, file)
            shutil.copy2(source_file, target_file)

    print(f"All .code-snippets files have been successfully moved to {target_dir}.")


# 清理临时目录
def clean_temp_directories():
    directories = [CONFIG["_download_dir"], CONFIG["_templates_dir"], CONFIG["_snippets_dir"]]
    for directory in directories:
        if os.path.exists(directory):
            shutil.rmtree(directory)
            print(f"Removed {directory} directory.")


# 主函数
def main():
    try:
        # 检查并安装必要的工具
        check_and_install_tools()
        # 从语雀下载文档
        download_yuque_doc()
        # 从 Markdown 文件提取代码到 templates 目录
        extract_code_to_templates()
        # 生成 code-snippets 文件
        generate_snippets()
        # 将 snippets 目录下的文件移动到目标目录
        move_snippets_to_target()
        # 清理临时目录
        clean_temp_directories()
        print("Successfully! The template library has been updated from Yuque to local snippets.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
```



> 更新: 2025-04-12 15:40:49  
> 原文: <https://www.yuque.com/capps/ctl/heb78p9y3xvyrpz9>