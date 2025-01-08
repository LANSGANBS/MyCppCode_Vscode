#include <stdio.h>
#include <stdlib.h>

typedef int ElemType;

typedef struct BiTNode
{
    ElemType data;
    struct BiTNode *left;
    struct BiTNode *right;
} BiTNode, *BiTree;

void InOrderTraverse(BiTree T)
{
    if (T != NULL)
    {
        InOrderTraverse(T->left);
        printf("%d ", T->data);
        InOrderTraverse(T->right);
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
    root->left->right = (BiTree)malloc(sizeof(BiTNode));
    root->left->right->data = 5;
    root->right = (BiTree)malloc(sizeof(BiTNode));
    root->right->data = 3;
    root->left->left->left = NULL;
    root->left->left->right = NULL;
    root->left->right->left = NULL;
    root->left->right->right = NULL;
    root->right->left = NULL;
    root->right->right = NULL;
    return root;
}

int main()
{
    BiTree tree = CreateSampleTree();
    printf("中序遍历结果: ");
    InOrderTraverse(tree);
    printf("\n");

    free(tree->left->right);
    free(tree->left->left);
    free(tree->left);
    free(tree->right);
    free(tree);

    return 0;
}
