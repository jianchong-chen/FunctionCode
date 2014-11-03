#ifndef _PLOYGON_FILL_H_
#define _PLOYGON_FILL_H_

#include "ai_defs.h"

#define MAX_PLOYGON_NUM 10
#define DEFAUlT_HEIGHT 1080
#define DEFAULT_WIDTH  1920

//基于扫描线的多边形填充算法
class CPloygonFill
{
public:
    CPloygonFill(l32 l32W = DEFAULT_WIDTH, l32 l32H = DEFAUlT_HEIGHT, TPolygon *ptPloygon = NULL, l32 l32Num = 0);
    ~CPloygonFill();

    const u8 *GetMask();
    l32 AddPloygon(TPolygon *ptPloygon, l32 l32Num = 1);
    void Reset();
    void SetParam(l32 l32H, l32 l32W);

protected:
    u8 *pu8Mask;
    l32 l32Width;
    l32 l32Height;

    TPolygon atPloygon[MAX_PLOYGON_NUM];
    l32 l32PloyNum;
};

#endif
