#include <stdio.h>
#include <stdlib.h>

typedef int ElemType;

typedef struct BiTNode
{
    ElemType data;
    struct BiTNode *left;
    struct BiTNode *right;
} BiTNode, *BiTree;

void PostOrderTraverse(BiTree T)
{
    if (T != NULL)
    {
        PostOrderTraverse(T->left);
        PostOrderTraverse(T->right);
        printf("%d ", T->data);
    }
}

void FreeTree(BiTree T)
{
    if (T != NULL)
    {
        FreeTree(T->left);
        FreeTree(T->right);
        free(T);
    }
}

BiTree CreateSampleTree()
{
    BiTree root = (BiTree)malloc(sizeof(BiTNode));
    root->data = 1;
    root->left = (BiTree)malloc(sizeof(BiTNode));
    root->left->data = 2;
    root->left->left = (BiTree)malloc(sizeof(BiTNode));
    root->left->left->data = 4;
    root->left->left->left = NULL;
    root->left->left->right = NULL;
    root->left->right = (BiTree)malloc(sizeof(BiTNode));
    root->left->right->data = 5;
    root->left->right->left = NULL;
    root->left->right->right = NULL;
    root->right = (BiTree)malloc(sizeof(BiTNode));
    root->right->data = 3;
    root->right->left = NULL;
    root->right->right = NULL;
    return root;
}

int main()
{
    BiTree tree = CreateSampleTree();
    printf("后序遍历结果: ");
    PostOrderTraverse(tree);
    printf("\n");

    FreeTree(tree);

    return 0;
}
