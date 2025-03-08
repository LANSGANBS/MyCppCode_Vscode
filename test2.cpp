#include <bits/stdc++.h>
using namespace std;

#define endl '\n'
#define all(x) begin(x), end(x)
#define mem(a, x) memset(a, x, sizeof(a))
#define int ll
typedef long long ll;

const int NMAX = 50000 + 5;

enum EventType { QUERY, MODIFY };

struct Query {
  int l, r, t, idx;
};

struct Modification {
  int pos;    // 修改位置
  int type;   // 1: 修改颜色, 2: 修改数字
  int prevC;  // 若是颜色修改，原来的颜色
  int nowC;
  int prevA;  // 若是数字修改，原来的数字
  int nowA;
};

int n, m, k;
int initA[NMAX], curA[NMAX];
int initC[NMAX], curC[NMAX];  // curC 存储离散化后的颜色
vector<int> compColors;

int totQuery = 0, totMod = 0;
vector<Query> queries;
vector<Modification> mods;

ll ansArr[NMAX];

// 优化：使用 vector 记录每个颜色（离散化后编号）的当前数字和
vector<ll> colSumArr;  // 下标范围 [1, totColors]

ll curAns = 0;

// 快速计算 f(x) = x^k，其中 k=1,2,3
inline ll f(ll x) {
  if (k == 1)
    return x;
  else if (k == 2)
    return x * x;
  else
    return x * x * x;
}

// 内联更新函数：移除位置 pos 的贡献
inline void removePos(int pos) {
  int color = curC[pos];  // 范围 1..totColors
  int val = curA[pos];
  ll oldVal = colSumArr[color];
  ll newVal = oldVal - val;
  curAns -= (f(oldVal) - f(newVal));
  colSumArr[color] = newVal;
}

// 内联更新函数：添加位置 pos 的贡献
inline void addPos(int pos) {
  int color = curC[pos];
  int val = curA[pos];
  ll oldVal = colSumArr[color];
  ll newVal = oldVal + val;
  curAns += (f(newVal) - f(oldVal));
  colSumArr[color] = newVal;
}

// 应用修改操作 modIdx（如果 pos 在当前区间内，先 remove 后 add）
void applyModification(int modIdx, int L, int R) {
  Modification &mod = mods[modIdx];
  int pos = mod.pos;
  if (L <= pos && pos <= R) removePos(pos);
  if (mod.type == 1) {
    // 颜色修改：curC[pos] 从 prevC 变为 nowC
    curC[pos] = mod.nowC;
  } else {
    curA[pos] = mod.nowA;
  }
  if (L <= pos && pos <= R) addPos(pos);
}

// 撤销修改操作 modIdx
void undoModification(int modIdx, int L, int R) {
  Modification &mod = mods[modIdx];
  int pos = mod.pos;
  if (L <= pos && pos <= R) removePos(pos);
  if (mod.type == 1) {
    curC[pos] = mod.prevC;
  } else {
    curA[pos] = mod.prevA;
  }
  if (L <= pos && pos <= R) addPos(pos);
}

void solve() {
  cin >> n >> m >> k;
  for (int i = 1; i <= n; i++) {
    cin >> initA[i];
    curA[i] = initA[i];
  }
  for (int i = 1; i <= n; i++) {
    cin >> initC[i];
    compColors.push_back(initC[i]);
    curC[i] = initC[i];  // 后续会进行离散化
  }

  ll lastans = 0;
  for (int i = 0; i < m; i++) {
    char op;
    cin >> op;
    if (op == 'Q') {
      int l, r;
      cin >> l >> r;
      l ^= lastans;
      r ^= lastans;
      if (l > r) swap(l, r);
      queries.push_back({l, r, (int)totMod, (int)totQuery});
      totQuery++;
    } else if (op == 'C') {
      int x, y;
      cin >> x >> y;
      x ^= lastans;
      y ^= lastans;
      compColors.push_back(y);
      mods.push_back({x, 1, curC[x], y, 0, 0});
      curC[x] = y;
      totMod++;
    } else {  // op == 'R'
      int x, y;
      cin >> x >> y;
      x ^= lastans;
      y ^= lastans;
      mods.push_back({x, 2, 0, 0, curA[x], y});
      curA[x] = y;
      totMod++;
    }
  }

  // 离散化：压缩所有颜色
  sort(compColors.begin(), compColors.end());
  compColors.erase(unique(compColors.begin(), compColors.end()),
                   compColors.end());
  int totColors = compColors.size();
  auto getColorId = [&](int c) -> int {
    return (int)(lower_bound(compColors.begin(), compColors.end(), c) -
                 compColors.begin()) +
           1;
  };
  for (int i = 1; i <= n; i++) {
    initC[i] = getColorId(initC[i]);
    curC[i] = initC[i];
  }
  for (auto &mod : mods) {
    if (mod.type == 1) {
      mod.prevC = getColorId(mod.prevC);
      mod.nowC = getColorId(mod.nowC);
    }
  }

  for (int i = 1; i <= n; i++) {
    curA[i] = initA[i];
    curC[i] = initC[i];
  }

  // 初始化 colSumArr，所有颜色贡献初始为 0
  colSumArr.assign(totColors + 1, 0);
  curAns = 0;

  // MO 排序：块大小取 n^(2/3)
  int blockSize = max((ll)1, (ll)pow(n, 2.0 / 3.0));
  sort(queries.begin(), queries.end(), [&](const Query &A, const Query &B) {
    int ablock = A.l / blockSize, bblock = B.l / blockSize;
    if (ablock != bblock) return ablock < bblock;
    int rblockA = A.r / blockSize, rblockB = B.r / blockSize;
    if (rblockA != rblockB) return A.r < B.r;
    return A.t < B.t;
  });

  // 恢复 curA, curC 为初始状态
  for (int i = 1; i <= n; i++) {
    curA[i] = initA[i];
    curC[i] = initC[i];
  }
  // colSumArr 已经初始化为 0，curAns = 0
  int L = 1, R = 0, curT = 0;
  for (auto &q : queries) {
    while (curT < q.t) {
      applyModification(curT, L, R);
      curT++;
    }
    while (curT > q.t) {
      curT--;
      undoModification(curT, L, R);
    }
    while (R < q.r) {
      R++;
      addPos(R);
    }
    while (R > q.r) {
      removePos(R);
      R--;
    }
    while (L < q.l) {
      removePos(L);
      L++;
    }
    while (L > q.l) {
      L--;
      addPos(L);
    }
    ansArr[q.idx] = curAns;
    lastans = curAns;
  }

  for (int i = 0; i < totQuery; i++) {
    cout << ansArr[i] << "\n";
  }
}

signed main() {
  // 重定向标准输入输出到指定文件
  freopen("C:\\Users\\70510\\Desktop\\05.in", "r", stdin);
  freopen("C:\\Users\\70510\\Desktop\\05.out", "w", stdout);

  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  int tt = 1;
  while (tt--) {
    solve();
  }
  return 0;
}
