#include <bits/stdc++.h>
using namespace std;
const int MAXN = 50050;
const int BLNB = 550;
const int COL = 1000050;
void read(int &x) {
  char ch;
  while (ch = getchar(), ch < '!');
  x = ch - 48;
  while (ch = getchar(), ch > '!') x = (x << 3) + (x << 1) + ch - 48;
}
int target[MAXN];
struct Change {
  int p, col, las;
} change[MAXN];
int nc, n, m, mp[COL], tot, D, cnt[BLNB][MAXN * 2], c[MAXN], blnm, spe[BLNB];
int CNT[MAXN * 2], ans[BLNB][BLNB], ima, tim[BLNB][BLNB], id[MAXN];
// 细节：我们不能直接在cnt[][]上做更改，所以需要记录一个临时的变化量数组CNT[]
// 变量解释：nc表示当前时间，mp[]和tot是离散化用的，D表示特征点步长，cnt[][]是预处理的莫队信息，id[]记录下标为i的特征点是第几个特征点，spe[]用于存储所有的特征点下标，ans[][]表示特征区间的答案，tim[][]记录答案的上一次更新时间，target[]表示离位置i最近的特征点坐标。
inline int getc(int sl, int sr, int p) {
  if (sl == 0 && sr == 0)
    return 0;
  else {
    if (c[sl] == p)
      return cnt[id[sr]][p] - cnt[id[sl]][p] + 1;
    else
      return cnt[id[sr]][p] - cnt[id[sl]][p];
  }  // 细节：端点特判一下
}
// 函数作用：读取区间[sl,sr]中的莫队数组信息。
inline void del(int pos, int sl, int sr) {
  if ((--CNT[c[pos]]) + getc(sl, sr, c[pos]) == 0) --ima;
}
inline void add(int pos, int sl, int sr) {
  if ((++CNT[c[pos]]) + getc(sl, sr, c[pos]) == 1) ++ima;
}
int main() {
  read(n);
  read(m);
  D = pow(n, 2.0 / 3);  // 带修莫队的块大小
  for (int i = 1; i <= n; ++i) {
    read(c[i]);
    if (!mp[c[i]])
      c[i] = mp[c[i]] = ++tot;
    else
      c[i] = mp[c[i]];
  }
  int tmp = 1;
  spe[blnm = 1] = 1;
  id[1] = 1;
  for (int i = 1; i <= n; ++i) {
    if (i - tmp == D) tmp = i, spe[++blnm] = i, id[i] = blnm;
    target[i] = tmp;
  }  // 预处理特征点以及每个点对应的离它最近的特征点
  int p = 1;
  for (int i = 1; i <= n; ++i) {
    ++CNT[c[i]];
    if (i == spe[p]) {
      for (int j = 1; j <= n; ++j) cnt[p][j] = CNT[j];
      ++p;
    }
  }  // 预处理莫队所需信息
  for (int i = 1; i <= blnm; ++i) {
    int p = i + 1;
    ima = 0;
    memset(CNT, 0, sizeof CNT);
    for (int j = spe[i]; j <= n; ++j) {
      if ((++CNT[c[j]]) == 1) ++ima;
      if (j == spe[p]) {
        ans[i][p] = ima;
        ++p;
      }
    }
  }  // 预处理特征区间答案
  memset(CNT, 0, sizeof CNT);
  while (m--) {
    char opt;
    int l, r;
    ima = 0;
    while (opt = getchar(), opt != 'Q' && opt != 'R');
    read(l);
    read(r);
    if (opt == 'R') {
      change[++nc].p = l;
      if (!mp[r])
        r = mp[r] = ++tot;
      else
        r = mp[r];
      change[nc].col = r;
      change[nc].las = c[l];
      int p = blnm;
      for (; spe[p] >= l; --p)
        --cnt[p][c[l]], ++cnt[p][r];  // 修改中间过程的信息
      c[l] = r;
    } else {
      int sl = target[l], sr = target[r];
      int SL = sl, SR = sr;
      // sl、sr表示所需特征区间的左右端点。
      if (sl == sr) {
        for (int i = l; i <= r; ++i)
          if (++CNT[c[i]] == 1) ++ima;
        printf("%d\n", ima);  // 细节：区间左右端点所属特征点相同，暴力计算
        for (int i = l; i <= r; ++i) --CNT[c[i]];
        // 临时数组还原
      } else {
        for (int t = nc; t > tim[id[sl]][id[sr]]; --t) {
          if (sl <= change[t].p && change[t].p <= sr) {
            if (++CNT[change[t].las] + getc(sl, sr, change[t].las) == 1)
              --ans[id[sl]][id[sr]];
            if (--CNT[change[t].col] + getc(sl, sr, change[t].col) == 0)
              ++ans[id[sl]][id[sr]];
            // 反向计算答案，即把原来带修莫队的东西反过来写，详情可以参考题解里普通带修莫队的修改方式做对照。
          }
        }
        for (int t = nc; t > tim[id[sl]][id[sr]]; --t)
          if (sl <= change[t].p && change[t].p <= sr) {
            --CNT[change[t].las];
            ++CNT[change[t].col];
          }
        // 临时数组还原
        tim[id[sl]][id[sr]] = nc;
        ima = ans[id[sl]][id[sr]];
        while (sl < l) del(sl++, SL, SR);
        while (sl > l) add(--sl, SL, SR);
        while (sr < r) add(++sr, SL, SR);
        while (sr > r) del(sr--, SL, SR);
        // 可爱的四句莫队
        printf("%d\n", ima);
        while (SL < l) add(SL++, 0, 0);
        while (SL > l) del(--SL, 0, 0);
        while (SR < r) del(++SR, 0, 0);
        while (SR > r) add(SR--, 0, 0);
        // 临时数组还原
      }
    }
  }
}