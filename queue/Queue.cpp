#include "Queue.h"
#include <stdlib.h>
#ifdef AVMalloc
#define MALLOC AVMalloc
#define FREE   AVFree
#else
#define MALLOC malloc
#define FREE   free
#endif

TQueue::TQueue(int l32Size)
{
    pData = (DataType *)MALLOC(l32Size * sizeof(DataType));
    l32MaxLen = l32Size;
    l32Start = 0;
    l32EleNum = 0;
}

TQueue::~TQueue()
{
    if(pData)
    {
        FREE(pData);
        pData = NULL;
    }
}

void TQueue::push(DataType tInfo)
{
    if(l32EleNum >= l32MaxLen) 
    {
        pop();//Ñ­»·¸²¸Ç
        //return; //²»¸²¸Ç
    }
    int l32Id = l32Start + l32EleNum;
    if(l32Id >= l32MaxLen) l32Id -= l32MaxLen;
    pData[l32Id] = tInfo;
    l32EleNum++;
};

void TQueue::pop()
{
    l32Start++;
    l32EleNum--;
    if(l32Start >= l32MaxLen) l32Start -= l32MaxLen;
    if(l32EleNum < 0) l32EleNum = 0;
};

DataType *TQueue::Front()
{
    return GetElem(0);
}

DataType *TQueue::Tail()
{
    return GetElem(l32EleNum - 1);
}

DataType *TQueue::GetElem(int l32Idx)
{
    if(l32Idx < 0 || l32Idx >= l32EleNum)
    {
        return NULL;
    }
    int l32InterID = l32Start + l32Idx;
    if(l32InterID >= l32MaxLen) 
        l32InterID -= l32MaxLen;
    return &pData[l32InterID];
}

DataType *TQueue::operator[](int l32Idx)
{
    if(l32Idx < 0 || l32Idx >= l32EleNum)
    {
        return NULL;
    }
    int l32InterID = l32Start + l32Idx;
    if(l32InterID >= l32MaxLen) 
        l32InterID -= l32MaxLen;
    return &pData[l32InterID];
}