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

i64 ceilDiv(i64 n, i64 m)
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

i64 floorDiv(i64 n, i64 m)
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

const int mod = 1e9 + 7;
constexpr int N = 2e5 + 7;
constexpr int M = 2e3 + 7;

typedef char ElemType;

typedef struct CSLinklist
{
    ElemType data;
    struct CSLinklist *firstChild;
    struct CSLinklist *nextSibling;
} CSLinklist, *CSTree;

CSTree createNode(ElemType data)
{
    CSTree newNode = (CSTree)malloc(sizeof(CSLinklist));
    newNode->data = data;
    newNode->firstChild = NULL;
    newNode->nextSibling = NULL;
    return newNode;
}

void printTree(CSTree tree, int level)
{
    if (tree)
    {
        for (int i = 0; i < level; i++)
        {
            printf("  ");
        }
        printf("%c\n", tree->data);
        printTree(tree->firstChild, level + 1);
        printTree(tree->nextSibling, level);
    }
}

void freeTree(CSTree tree)
{
    if (tree)
    {
        freeTree(tree->firstChild);
        freeTree(tree->nextSibling);
        free(tree);
    }
}

void deleteNode(CSTree *tree, ElemType target)
{
    if (!(*tree))
        return;

    if ((*tree)->data == target)
    {
        freeTree(*tree);
        *tree = NULL;
        return;
    }
    if ((*tree)->firstChild && (*tree)->firstChild->data == target)
    {
        freeTree((*tree)->firstChild);
        (*tree)->firstChild = NULL;
        return;
    }
    CSTree sibling = (*tree)->firstChild;
    while (sibling && sibling->nextSibling)
    {
        if (sibling->nextSibling->data == target)
        {
            freeTree(sibling->nextSibling);
            sibling->nextSibling = NULL;
            return;
        }
        sibling = sibling->nextSibling;
    }
    deleteNode(&(*tree)->firstChild, target);
    deleteNode(&(*tree)->nextSibling, target);
}

void solve()
{
    CSTree root = createNode('T');
    root->firstChild = createNode('B');
    root->firstChild->nextSibling = createNode('D');
    root->firstChild->firstChild = createNode('E');
    root->firstChild->firstChild->nextSibling = createNode('F');
    root->firstChild->nextSibling->firstChild = createNode('G');
    printf("原始家谱树结构:\n");
    printTree(root, 0);
    printf("\n删除节点 'B' 后的家谱树结构:\n");
    deleteNode(&root, 'B');
    printTree(root, 0);
    freeTree(root);
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