#ifndef __QUEUE_H__
#define __QUEUE_H__
#define DEFAULT_LEN 256
typedef int DataType;

#include <stdio.h>

typedef struct TQueue
{
    DataType *pData;
    int l32MaxLen;
    int l32Start;
    int l32EleNum;

    TQueue(int l32Size = DEFAULT_LEN);
    ~TQueue();
    void push(DataType tInfo);
    void pop();
    DataType *Front();
    DataType *Tail();
    DataType *GetElem(int l32Idx); 
    DataType *operator[](int l32Idx);
};

#endif
