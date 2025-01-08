#include <bits/stdc++.h>
using namespace std;
using i64 = long long;
const int N = 2.01e5;
const int M = 1.01e3;
const int mod = 998244353;
#define endl '\n'
const int inf = 0x7fffffff;

int n;
int grid[12][12][12];
vector<vector<tuple<int, int, int>>> sorted_choices;
int min_sum = 0x7fffffff;
int min_possible[12];
unordered_map<int, int> memo;
chrono::steady_clock::time_point start_time;

#pragma GCC target("avx")
#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-fwhole-program")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-fstrict-overflow")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-skip-blocks")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("-funsafe-loop-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")

void solve()
{
    cin >> n;
    for (int x = 0; x < n; ++x)
    {
        for (int y = 0; y < n; ++y)
        {
            for (int z = 0; z < n; ++z)
            {
                cin >> grid[x][y][z];
            }
        }
    }
    sorted_choices.assign(n, vector<tuple<int, int, int>>());
    for (int x = 0; x < n; ++x)
    {
        for (int y = 0; y < n; ++y)
        {
            for (int z = 0; z < n; ++z)
            {
                sorted_choices[x].emplace_back(y, z, grid[x][y][z]);
            }
        }
        sort(sorted_choices[x].begin(), sorted_choices[x].end(), [&](const tuple<int, int, int> &a, const tuple<int, int, int> &b) -> bool
             { return get<2>(a) < get<2>(b); });
        min_possible[x] = get<2>(sorted_choices[x][0]);
    }
    min_sum = inf;
    start_time = chrono::steady_clock::now();
    function<void(int, int, int, int)> dfs = [&](int x, int y_mask, int z_mask, int current_sum) -> void
    {
        auto current_time = chrono::steady_clock::now();
        auto elapsed_time = chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count();
        if (elapsed_time > 3000)
        {
            return;
        }
        if (x == n)
        {
            if (current_sum < min_sum)
            {
                min_sum = current_sum;
            }
            return;
        }
        int state = (x << 24) | (y_mask << 12) | z_mask;
        if (memo.count(state) && memo[state] <= current_sum)
        {
            return;
        }
        memo[state] = current_sum;
        if (current_sum + min_possible[x] >= min_sum)
        {
            return;
        }
        for (auto &[y, z, val] : sorted_choices[x])
        {
            if (!(y_mask & (1 << y)) && !(z_mask & (1 << z)))
            {
                int new_sum = current_sum + val;
                if (new_sum >= min_sum)
                {
                    continue;
                }
                dfs(x + 1, y_mask | (1 << y), z_mask | (1 << z), new_sum);
            }
        }
    };
    dfs(0, 0, 0, 0);
    if (min_sum < inf)
    {
        cout << min_sum << endl;
    }
    else
    {
        cout << -1 << endl;
    }
}

signed main()
{
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    int T = 1;
    // cin >> T;
    while (T--)
        solve();
}