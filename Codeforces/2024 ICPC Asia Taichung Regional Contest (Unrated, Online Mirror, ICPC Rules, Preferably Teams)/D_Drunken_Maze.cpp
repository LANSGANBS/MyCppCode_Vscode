#include <bits/stdc++.h>
using namespace std;
#define int long long
using i64 = long long;
const int N_MAX = 2.01e5;
const int M_MAX = 1.01e3;
const int mod = 998244353;
#define endl '\n'

struct State
{
    int x, y;
    int last_dir;
    int count_steps;
};

void solve()
{
    int n, m;
    cin >> n >> m;
    vector<string> maze(n);
    int start_x = -1, start_y = -1;
    int end_x = -1, end_y = -1;
    for (int i = 0; i < n; ++i)
    {
        cin >> maze[i];
        for (int j = 0; j < m; ++j)
        {
            if (maze[i][j] == 'S')
            {
                start_x = i;
                start_y = j;
            }
            if (maze[i][j] == 'T')
            {
                end_x = i;
                end_y = j;
            }
        }
    }

    if (start_x == -1 || end_x == -1)
    {
        cout << -1 << endl;
        return;
    }

    vector<pair<int, int>> directions = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}};

    int total_cells = n * m;
    vector<bool> visited(total_cells * 12, false);

    auto get_cell_id = [&](int x, int y) -> int
    {
        return x * m + y;
    };

    queue<pair<State, int>> q;

    State start_state = {start_x, start_y, -1, 0};
    q.push({start_state, 0});

    bool found = false;
    int result_steps = -1;

    while (!q.empty())
    {
        auto current = q.front();
        q.pop();
        State s = current.first;
        int steps = current.second;

        if (s.x == end_x && s.y == end_y)
        {
            result_steps = steps;
            found = true;
            break;
        }

        for (int d = 0; d < 4; ++d)
        {
            int nx = s.x + directions[d].first;
            int ny = s.y + directions[d].second;
            if (nx < 0 || nx >= n || ny < 0 || ny >= m)
                continue;
            if (maze[nx][ny] == '#')
                continue;

            int new_dir = d;
            int new_count_steps;

            if (s.last_dir == d)
            {
                if (s.count_steps == 3)
                {
                    continue;
                }
                new_count_steps = s.count_steps + 1;
            }
            else
            {
                new_dir = d;
                new_count_steps = 1;
            }
            int cell_id = get_cell_id(nx, ny);
            int step_index = new_count_steps - 1;
            int state_idx;
            if (s.last_dir == -1)
            {
                state_idx = d * 3 + 0;
            }
            else
            {
                state_idx = new_dir * 3 + step_index;
            }

            int unique_id = cell_id * 12 + state_idx;

            if (!visited[unique_id])
            {
                visited[unique_id] = true;
                State new_state = {nx, ny, new_dir, new_count_steps};
                q.push({new_state, steps + 1});
            }
        }
    }

    cout << result_steps << endl;
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