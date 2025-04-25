下面是“缓存系统”问题的完整实现。思路：

对每个题目 (i) 预处理前缀：

$(\text{ws}[t]=\sum_{j=1}^t S_{ij})$；

$(\text{vs}[t]=\sum_{j=1}^t A_{ij})$。

总读次数 $(\text{total}=\sum_{i,j}A_{ij})$。

做多选背包：每道题目视作一组，组内可选前缀长度 $(t\in[0,M])$（$(t=0)$即不选任何数据），权重 $(\text{ws}[t])$，收益 $(\text{vs}[t])$。
背包容量为 $(X)$，求最大收益。答案就是 $(\text{total}-\text{maxSaved})$。