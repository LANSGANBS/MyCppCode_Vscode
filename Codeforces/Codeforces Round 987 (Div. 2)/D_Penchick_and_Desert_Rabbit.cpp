#include <bits/stdc++.h>

using namespace std;

#define endl "\n"
#define i64 long long
const int N = 1e6 + 10;

struct node
{
    int mx1, mx2;
    node operator+(const node &b)
    {
        node c = b;
        if (mx1 > c.mx1)
            c.mx2 = c.mx1, c.mx1 = mx1;
        else if (mx1 > c.mx2)
            c.mx2 = mx1;
        if (mx1 > c.mx2)
            c.mx2 = mx2;
        return c;
    }
    void op(int x)
    {
        mx1 = x;
        mx2 = 0;
    }
};
struct
{
    int l, r;
    node c;
    void op(int x)
    {
        c.op(x);
    }
} tr[N << 2];

void build(int u, int l, int r)
{
    tr[u].l = l;
    tr[u].r = r;
    tr[u].c.op(0);
    if (l == r)
        return;
    int mid = l + r >> 1;
    build(u << 1 | 1, mid + 1, r);
    build(u << 1, l, mid);
}
void push_up(int u)
{
    tr[u].c = tr[u << 1].c + tr[u << 1 | 1].c;
}
void modify(int u, int l, int r, int d)
{
    if (tr[u].l >= l && tr[u].r <= r)
    {
        tr[u].op(d);
        return;
    }
    else
    {
        int mid = tr[u].l + tr[u].r >> 1;
        if (l <= mid)
            modify(u << 1, l, r, d);
        if (r > mid)
            modify(u << 1 | 1, l, r, d);
        push_up(u);
    }
}
node ask(int u, int l, int r)
{
    if (tr[u].l >= l && tr[u].r <= r)
        return tr[u].c;
    int mid = tr[u].l + tr[u].r >> 1;
    if (r <= mid)
        return ask(u << 1, l, r);
    if (l > mid)
        return ask(u << 1 | 1, l, r);
    return ask(u << 1, l, r) + ask(u << 1 | 1, l, r);
}

void solve()
{
    int n;
    cin >> n;
    build(1, 1, n);
    vector<i64> pre(n), a(n);
    for (int i = 0; i < n; i++)
        cin >> a[i];
    pre[0] = a[0];
    for (int i = 1; i < n; i++)
        pre[i] = max(pre[i - 1], a[i]);
    vector<int> ans(n);
    for (int i = n - 1; ~i; i--)
    {
        ans[i] = pre[i];
        for (int k = 0; k < 3 && ans[i] > 1; k++)
            ans[i] = max(ans[i], ask(1, 1, ans[i] - 1).mx1);
        modify(1, a[i], a[i], ans[i]);
    }
    for (int i = 0; i < n; i++)
        cout << ans[i] << " \n"[i == n - 1];
}

signed main()
{
    cin.tie(0);
    cout.tie(0);
    ios::sync_with_stdio(0);
    int t;
    cin >> t;
    while (t--)
        solve();
}