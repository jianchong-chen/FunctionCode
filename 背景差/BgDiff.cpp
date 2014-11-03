#include "BgDiff.h"
#include <stdlib.h>
#include <string.h>

l32 BgDiffOpen(void **ppvFrameDiffHandle, l32 l32Height, l32 l32Width, l32 l32MinGrayDiff)
{
    TBgDiffHandle *ptHandle = NULL;
    ptHandle = (TBgDiffHandle *)malloc(sizeof(TBgDiffHandle));
    if(NULL == ptHandle)
    {
        goto ErrLabel;
    }
    memset(ptHandle, 0, sizeof(TBgDiffHandle));

    ptHandle->l32Width = l32Width;
    ptHandle->l32Height = l32Height;
    ptHandle->l32MinGrayDiff = l32MinGrayDiff;
    ptHandle->l32FrameNum = 0;
    
    ptHandle->pf32Buf = (f32 *)malloc(l32Height * l32Width * sizeof(f32));
    ptHandle->pf32Bg = (f32 *)malloc( l32Height * l32Width * sizeof(f32));
    if(NULL == ptHandle->pf32Buf || NULL == ptHandle->pf32Bg)
    {
        goto ErrLabel;
    }

    *ppvFrameDiffHandle = ptHandle;
    return 0;
ErrLabel:
    BgDiffClose(ptHandle);
    return 1;
}

l32 BgDiffProcess(void *pvFrameDiffHandle, u8 *pu8Img, u8 *pu8Fg)
{
    TBgDiffHandle *ptHandle = NULL;
    f32 *pf32Bg = NULL;
    f32 *pf32Diff = NULL;
    l32 l32Index;
    l32 l32Height, l32Width;  
    f32 f32Mean = 0.0f;
    f32 f32Th;
    f32 f32Rate = 0;
    f32 f32FgRate;
    f32 f32BgRate;

    if(NULL == pvFrameDiffHandle || NULL == pu8Img || NULL == pu8Fg)
    {
        goto End;
    }
    ptHandle = (TBgDiffHandle *)pvFrameDiffHandle;

    l32Height = ptHandle->l32Height;
    l32Width = ptHandle->l32Width;
    pf32Bg = ptHandle->pf32Bg;
    pf32Diff = ptHandle->pf32Buf;

    if(ptHandle->l32FrameNum == 0)
    {
        for(l32Index = 0; l32Index < l32Height * l32Width; l32Index++)
        {
            pf32Bg[l32Index] = (f32)pu8Img[l32Index];
        }
        memset(pu8Fg, DIFF_BG, l32Height * l32Width);
        goto End;
    }


    f32Mean = 0.0f;
    for(l32Index = 0; l32Index < l32Height * l32Width; l32Index++)
    {
        pf32Diff[l32Index] = ABS((f32)pu8Img[l32Index] - pf32Bg[l32Index]);
        f32Mean += pf32Diff[l32Index];
    }

    f32Th = 2.0f * f32Mean / (l32Height * l32Width);
    f32Th = MAX(f32Th, ptHandle->l32MinGrayDiff);

    f32BgRate = 1.0f / MIN(ptHandle->l32FrameNum , 500);
    f32FgRate = f32BgRate / 3;
    for(l32Index = 0; l32Index < l32Height * l32Width; l32Index++)
    {
        if(pf32Diff[l32Index] > f32Th) 
        {
            pu8Fg[l32Index] = DIFF_FG;
            f32Rate = f32FgRate;
        }
        else
        {
            pu8Fg[l32Index] = DIFF_BG;
            f32Rate = f32BgRate;
        }
        pf32Bg[l32Index] = pf32Bg[l32Index] * (1 - f32Rate) + pu8Img[l32Index] * f32Rate;
    }

End:
    ptHandle->l32FrameNum++;
    return 0;
}

l32 BgDiffClose(void *pvFrameDiffHandle)
{
    TBgDiffHandle *ptHandle = NULL;
    if(NULL == pvFrameDiffHandle)
    {
        goto End;
    }
    ptHandle = (TBgDiffHandle *)pvFrameDiffHandle;
    if(NULL != ptHandle->pf32Bg)
    {
        free(ptHandle->pf32Bg);
    }
    if(NULL != ptHandle->pf32Buf)
    {
        free(ptHandle->pf32Buf);
    }
    free(pvFrameDiffHandle);

End:
    return 0;

}
