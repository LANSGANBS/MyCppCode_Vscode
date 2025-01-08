#include <stdio.h>
#include <stdlib.h>

typedef int ElemType;

typedef struct BiTNode
{
    ElemType data;
    struct BiTNode *lchild;
    struct BiTNode *rchild;
} BiTNode, *BiTree;

typedef enum
{
    OK,
    ERROR
} Status;

Status PrintElement(ElemType e)
{
    printf("%d ", e);
    return OK;
}

Status PreOrderTraverse(BiTree T, Status (*Visit)(ElemType))
{
    if (T)
    {
        if (Visit(T->data) == ERROR)
            return ERROR;
        if (PreOrderTraverse(T->lchild, Visit) == ERROR)
            return ERROR;
        if (PreOrderTraverse(T->rchild, Visit) == ERROR)
            return ERROR;
    }
    return OK;
}

BiTree CreateSampleTree()
{
    BiTree root = (BiTree)malloc(sizeof(BiTNode));
    root->data = 1;
    root->lchild = (BiTree)malloc(sizeof(BiTNode));
    root->lchild->data = 2;
    root->lchild->lchild = NULL;
    root->lchild->rchild = NULL;
    root->rchild = (BiTree)malloc(sizeof(BiTNode));
    root->rchild->data = 3;
    root->rchild->lchild = NULL;
    root->rchild->rchild = NULL;
    return root;
}

int main()
{
    BiTree tree = CreateSampleTree();
    printf("前序遍历结果: ");
    PreOrderTraverse(tree, PrintElement);
    printf("\n");

    free(tree->lchild);
    free(tree->rchild);
    free(tree);

    return 0;
}
