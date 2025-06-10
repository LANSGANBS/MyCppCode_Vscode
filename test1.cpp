#include <bits/stdc++.h>
using namespace std;

char posChar(int idx) { return idx == 0 ? 'L' : (idx == 1 ? 'M' : 'R'); }

struct State {
  array<string, 3> comb;
  string key() const { return comb[0] + "|" + comb[1] + "|" + comb[2]; }
};

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  int T;
  cin >> T;
  while (T--) {
    int fig[3];
    for (int i = 0; i < 3; i++) {
      cin >> fig[i];
    }
    auto makeTarget = [&](int forbidden) -> string {
      vector<char> v;
      for (char d : {'0', '3', '4'}) {
        if (d - '0' != forbidden) v.push_back(d);
      }
      sort(v.begin(), v.end());
      return string(v.begin(), v.end());
    };
    array<string, 3> target;
    target[0] = makeTarget(fig[0]);
    target[1] = makeTarget(fig[1]);
    target[2] = makeTarget(fig[2]);
    State init;
    for (int i = 0; i < 3; i++) {
      string s;
      cin >> s;
      sort(s.begin(), s.end());
      init.comb[i] = s;
    }
    using Path = vector<string>;
    queue<pair<State, Path>> qu;
    unordered_set<string> vis;
    qu.push({init, {}});
    vis.insert(init.key());

    bool found = false;
    Path ans;

    while (!qu.empty() && !found) {
      auto cur = qu.front();
      qu.pop();
      State s = cur.first;
      Path path = cur.second;
      if (!path.empty() && s.comb[0] == target[0] && s.comb[1] == target[1] &&
          s.comb[2] == target[2]) {
        ans = path;
        found = true;
        break;
      }
      vector<tuple<string, State, Path>> nextStates;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          if (i == j) continue;
          for (char x : s.comb[i]) {
            for (char y : s.comb[j]) {
              if (x == y) continue;
              State ns = s;
              {
                string temp = ns.comb[i];
                for (auto &ch : temp) {
                  if (ch == x) {
                    ch = y;
                    break;
                  }
                }
                sort(temp.begin(), temp.end());
                ns.comb[i] = temp;
              }
              {
                string temp = ns.comb[j];
                for (auto &ch : temp) {
                  if (ch == y) {
                    ch = x;
                    break;
                  }
                }
                sort(temp.begin(), temp.end());
                ns.comb[j] = temp;
              }
              string op;
              op.push_back(x);
              op += "->";
              op.push_back(posChar(i));
              op += " ";
              op.push_back(y);
              op += "->";
              op.push_back(posChar(j));

              Path np = path;
              np.push_back(op);

              nextStates.push_back({op, ns, np});
            }
          }
        }
      }
      sort(nextStates.begin(), nextStates.end(),
           [](auto &a, auto &b) { return get<0>(a) < get<0>(b); });

      for (auto &entry : nextStates) {
        State ns = get<1>(entry);
        string key = ns.key();
        if (vis.count(key)) continue;
        vis.insert(key);
        qu.push({ns, get<2>(entry)});
      }
    }
    for (auto &s : ans) cout << s << "\n";
    cout << "---------" << endl;
  }
  return 0;
}