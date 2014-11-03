#include "PolygonFill.h"
#include <string.h>
#include <vector>
#include <list>
#include <algorithm>
using std::vector;
using std::list;

typedef struct
{
    l32 l32yMax;
    f32 f32x;
    f32 f32Dx;
}TEdge;

static void GetYMinMax(TPolygon *ptPloyGon, l32 *pl32Min, l32 *pl32Max);
static void InitEdge(TPolygon *ptPloyGon, vector< list<TEdge> > &vlAET);
static void UpdateAET(list<TEdge> &listAet, l32 l32Y);
static void AddToAET(list<TEdge> &listAet, list<TEdge> &NewAet);
static void FillLine(list<TEdge> &listAet, l32 l32Y, u8 *pu8Mask, l32 l32Stride);
static void PloygonFill(u8 *pu8Mask, l32 l32Height, l32 l32Width, TPolygon *ptPloyGon);


CPloygonFill::CPloygonFill(l32 l32W /* = DEFAULT_WIDTH */, l32 l32H /* = DEFAUlT_HEIGHT */, TPolygon *ptPloygon /* = NULL */, l32 l32Num /* = 0 */)
{
    l32Height = l32H;
    l32Width = l32W;
    l32PloyNum = 0;

    pu8Mask = new u8[l32H * l32W];
    memset(pu8Mask, 0, l32H * l32W * sizeof(u8));


    AddPloygon(ptPloygon, l32Num);
}

CPloygonFill::~CPloygonFill()
{
    l32PloyNum = 0;
    if(pu8Mask) 
    {
        delete []pu8Mask;
        pu8Mask = NULL;
    }
}

const u8 *CPloygonFill::GetMask()
{
    return pu8Mask;
}

l32 CPloygonFill::AddPloygon(TPolygon *ptPloygon, l32 l32Num /* = 1 */)
{
    if(l32PloyNum + l32Num > MAX_PLOYGON_NUM)
    {
        return 0;
    }
    if(l32Num == 0 || ptPloygon == NULL)
    {
        return 0;
    }

    for(l32 l32Index = 0; l32Index < l32Num; l32Index++)
    {
        atPloygon[l32Index + l32PloyNum] = ptPloygon[l32Index];

        PloygonFill(pu8Mask, l32Height, l32Width, ptPloygon + l32Index);
    }

    l32PloyNum += l32Num;
    
    return 1;
}

void CPloygonFill::Reset()
{
    l32PloyNum = 0;
    memset(pu8Mask, 0, l32Height * l32Width * sizeof(u8));
}

void CPloygonFill::SetParam(l32 l32H, l32 l32W)
{
    if(l32H * l32W > l32Height * l32Width)
    {
        if(pu8Mask) delete[]pu8Mask;
        pu8Mask = new u8[l32H * l32W];
        memset(pu8Mask, 0, l32H * l32W * sizeof(u8));
    }
    l32Height = l32H;
    l32Width = l32W;
}

static void PloygonFill(u8 *pu8Mask, l32 l32Height, l32 l32Width, TPolygon *ptPloyGon)
{
    vector< list<TEdge> > vlAET(l32Height);

    l32 l32yMin, l32yMax;
    GetYMinMax(ptPloyGon, &l32yMin, &l32yMax);

    InitEdge(ptPloyGon, vlAET);

    //开始扫描填充
    list<TEdge> listAet;
    for(l32 l32Index = l32yMin; l32Index <= l32yMax; l32Index++)
    {
        UpdateAET(listAet, l32Index);
        AddToAET(listAet, vlAET[l32Index]);
        FillLine(listAet, l32Index, pu8Mask, l32Width);
    }

    //填充水平线
    for(l32 l32Index = 0; l32Index < ptPloyGon->l32PointNum - 1; l32Index++)
    {
        l32 l32NxtId = l32Index + 1;
        if(l32NxtId > ptPloyGon->l32PointNum) l32NxtId -= ptPloyGon->l32PointNum;
        if(ptPloyGon->atPoint[l32Index].l32Y != ptPloyGon->atPoint[l32NxtId].l32Y)
        {
            continue;
        }
        l32 l32Y  = ptPloyGon->atPoint[l32Index].l32Y;
        l32 l32X1 = MIN(ptPloyGon->atPoint[l32Index].l32X, ptPloyGon->atPoint[l32Index + 1].l32X);
        l32 l32X2 = MAX(ptPloyGon->atPoint[l32Index].l32X, ptPloyGon->atPoint[l32Index + 1].l32X);
        memset(pu8Mask + l32Y * l32Width + l32X1, 1, l32X2 - l32X1 + 1);
    }
}

static void InitEdge(TPolygon *ptPloyGon, vector< list<TEdge> > &vlAET)
{
    //判断最后一个点是否与第一个点相同
    l32 l32PtNum = ptPloyGon->l32PointNum;
    TPoint *ptPoint = ptPloyGon->atPoint;
    if(ptPoint[l32PtNum - 1].l32X == ptPoint[0].l32X && ptPoint[l32PtNum - 1].l32Y == ptPoint[0].l32Y)
    {
        l32PtNum--;
    }

    //初始化活性边表
    for(l32 l32Index = 0; l32Index < l32PtNum; l32Index++)
    {
        l32 l32PreId = l32Index - 1;
        l32 l32NxtId = l32Index + 2;
        TPoint &tPCur1 = ptPoint[l32Index];
        TPoint &tPCur2 = ptPoint[l32Index + 1];
        if(l32PreId < 0) l32PreId += l32PtNum;
        if(l32NxtId >= l32PtNum) l32NxtId -= l32PtNum;

        TPoint &tPPre = ptPoint[l32PreId];
        TPoint &tPNxt = ptPoint[l32NxtId];
        
        //水平线
        if(tPCur1.l32Y == tPCur2.l32Y)
        {
            continue;
        }

        TEdge tEdge;
        l32 l32yMin = 0;
        
        tEdge.f32Dx = (tPCur1.l32X - tPCur2.l32X) * 1.0f / (tPCur1.l32Y - tPCur2.l32Y);
        if(tPCur1.l32Y > tPCur2.l32Y)
        {
            tEdge.f32x = (f32)tPCur2.l32X;
            tEdge.l32yMax = tPCur1.l32Y;
            l32yMin = tPCur2.l32Y;
            if(tPPre.l32Y >= tPCur1.l32Y)//左右顶点，只计算一次
            {
                tEdge.l32yMax--;
            }
        }
        else
        {
            tEdge.f32x = (f32)tPCur1.l32X;
            tEdge.l32yMax = tPCur2.l32Y;
            l32yMin = tPCur1.l32Y;
            if(tPNxt.l32Y >= tPCur2.l32Y)//左右顶点，只计算一次
            {
                tEdge.l32yMax--;
            }
        }
        vlAET[l32yMin].push_back(tEdge);
    }
}

static void GetYMinMax(TPolygon *ptPloyGon, l32 *pl32Min, l32 *pl32Max)
{
    l32 l32Min, l32Max;
    l32Min = l32Max = ptPloyGon->atPoint[0].l32Y;
    for(l32 l32Index = 1; l32Index < ptPloyGon->l32PointNum; l32Index++)
    {
        if(ptPloyGon->atPoint[l32Index].l32Y < l32Min)
        {
            l32Min = ptPloyGon->atPoint[l32Index].l32Y;
        }

        if(ptPloyGon->atPoint[l32Index].l32Y > l32Max)
        {
            l32Max = ptPloyGon->atPoint[l32Index].l32Y;
        }
    }

    *pl32Max = l32Max;
    *pl32Min = l32Min;
}
static void UpdateAET(list<TEdge> &listAet, l32 l32Y)
{
    list<TEdge>::iterator it;
    for(it = listAet.begin(); it != listAet.end();)
    {
        if(it->l32yMax < l32Y) 
        {
            it = listAet.erase(it);
        }
        else 
        {
            it->f32x += it->f32Dx;
            it++;
        }
    }
}

static void AddToAET(list<TEdge> &listAet, list<TEdge> &NewAet)
{
    list<TEdge>::iterator it;
    for(it = NewAet.begin(); it != NewAet.end(); it++)
    {
        listAet.push_back(*it);
    }
}

static void FillLine(list<TEdge> &listAet, l32 l32Y, u8 *pu8Mask, l32 l32Stride)
{
    vector<l32> vl32x;
    for(list<TEdge>::iterator it = listAet.begin(); it != listAet.end(); it++)
    {
        vl32x.push_back((l32)(it->f32x + 0.5));
    }
    sort(vl32x.begin(), vl32x.end());
    pu8Mask = pu8Mask + l32Y * l32Stride;
    for(l32 l32Index = 0; l32Index < (l32)vl32x.size() - 1; l32Index += 2)
    {
        l32 l32x1 = vl32x[l32Index];
        l32 l32x2 = vl32x[l32Index + 1];       
        memset(pu8Mask + l32x1, 1, l32x2 - l32x1 + 1);        
    }
}