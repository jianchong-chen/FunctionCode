#include "FrameDiff.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void CalcStd(u8 *pu8Gray, f32 *pf32Contrast, l32 * pl32IImg, l32 *pl32IImg2, l32 l32Width, l32 l32Height, l32 l32BlockW, l32 l32BlockH);

l32 FrameDiffOpen(void **ppvFrameDiffHandle, l32 l32Height, l32 l32Width, l32 l32DiffFrameNum, l32 l32MinGrayDiff, l32 l32MinStdDiff)
{
    TFrameDiffHandle *ptHandle = NULL;
    l32 l32ImgBufNum, l32Index;
    
    ptHandle = (TFrameDiffHandle *)malloc(sizeof(TFrameDiffHandle));
    if(NULL == ptHandle)
    {
        goto ErrLabel;
    }
    memset(ptHandle, 0, sizeof(TFrameDiffHandle));

    ptHandle->l32Width = l32Width;
    ptHandle->l32Height = l32Height;
    ptHandle->l32DiffFrameNum = l32DiffFrameNum;
    ptHandle->l32MinGrayDiff = l32MinGrayDiff;
    ptHandle->l32MinStdDiff = l32MinStdDiff;
    ptHandle->l32FrameNum = 0;
    l32ImgBufNum = 2 * ptHandle->l32DiffFrameNum + 1;

    ptHandle->pu8ImgBuf = (u8 *)malloc(l32ImgBufNum * l32Height * l32Width * sizeof(u8));
    ptHandle->pf32StdBuf = (f32 *)malloc(l32ImgBufNum * l32Height * l32Width * sizeof(f32));
    ptHandle->ppu8Img = (u8 **)malloc(l32ImgBufNum * sizeof(u8 *));
    ptHandle->ppf32Std = (f32 **)malloc(l32ImgBufNum * sizeof(f32 *));
    if(NULL == ptHandle->pu8ImgBuf || NULL == ptHandle->pf32StdBuf || NULL == ptHandle->ppu8Img || NULL == ptHandle->ppf32Std)
    {
        goto ErrLabel;
    }

    for(l32Index = 0; l32Index < l32ImgBufNum; l32Index++)
    {
        (ptHandle->ppu8Img)[l32Index] = ptHandle->pu8ImgBuf + l32Index * l32Height * l32Width;
        (ptHandle->ppf32Std)[l32Index] = ptHandle->pf32StdBuf + l32Index * l32Height * l32Width;
    }

    ptHandle->pvBuffer = malloc(2 * (l32Height + 1) * (l32Width + 1) * sizeof(f32));
    if(NULL == ptHandle->pvBuffer)
    {
        goto ErrLabel;
    }

    *ppvFrameDiffHandle = ptHandle;
    return 0;
ErrLabel:
    FrameDiffClose(ptHandle);
    return 1;
}

l32 FrameDiffProcess(void *pvFrameDiffHandle, u8 *pu8Img, u8 *pu8Fg)
{
    TFrameDiffHandle *ptHandle = NULL;
    l32 l32PreIdx, l32CurIdx, l32NxtIdx, l32ImgBufNum, l32Index;
    l32 l32Height, l32Width;
    u8 *pu8Pre = NULL;
    u8 *pu8Cur = NULL;
    u8 *pu8Nxt = NULL;
    f32 *pf32Pre = NULL;
    f32 *pf32Cur = NULL;
    f32 *pf32Nxt = NULL;
    f32 *pf32Buf1 = NULL;
    f32 *pf32Buf2 = NULL;
    f32 *pf32GD = NULL;
    f32 *pf32SD = NULL;
    f32 f32GMean = 0.0f;
    f32 f32SMean = 0.0f;
    f32 f32GTh, f32STh;

    if(NULL == pvFrameDiffHandle || NULL == pu8Img || NULL == pu8Fg)
    {
        goto End;
    }
    ptHandle = (TFrameDiffHandle *)pvFrameDiffHandle;

    l32ImgBufNum = 2 * ptHandle->l32DiffFrameNum + 1;

    l32Height = ptHandle->l32Height;
    l32Width = ptHandle->l32Width;

    l32PreIdx = ptHandle->l32FrameNum % l32ImgBufNum;
    l32CurIdx = (ptHandle->l32FrameNum + ptHandle->l32DiffFrameNum) % l32ImgBufNum;
    l32NxtIdx = (ptHandle->l32FrameNum + 2 * ptHandle->l32DiffFrameNum) % l32ImgBufNum;

    pu8Pre = ptHandle->ppu8Img[l32PreIdx];
    pu8Cur = ptHandle->ppu8Img[l32CurIdx];
    pu8Nxt = ptHandle->ppu8Img[l32NxtIdx];

    pf32Pre = ptHandle->ppf32Std[l32PreIdx];
    pf32Cur = ptHandle->ppf32Std[l32CurIdx];
    pf32Nxt = ptHandle->ppf32Std[l32NxtIdx];

    pf32Buf1 = (f32 *)ptHandle->pvBuffer;
    pf32Buf2 = pf32Buf1 + (l32Height + 1) * (l32Width + 1);

    memcpy(pu8Nxt, pu8Img, l32Height * l32Width * sizeof(u8));
    CalcStd(pu8Nxt, pf32Nxt, (l32 *)pf32Buf1, (l32 *)pf32Buf2, l32Width, l32Height, 5, 5);

    if(ptHandle->l32FrameNum < l32ImgBufNum - 1)
    {
        memset(pu8Fg, FRAME_DIFF_BG, ptHandle->l32Width * ptHandle->l32Height);
        goto End;
    }
    pf32GD = pf32Buf1;
    pf32SD = pf32Buf2;
    f32GMean = 0.0f;
    f32SMean = 0.0f;
    for(l32Index = 0; l32Index < l32Height * l32Width; l32Index++)
    {
        l32 l32GD1 = (l32)pu8Cur[l32Index] - (l32)pu8Pre[l32Index];
        l32 l32GD2 = (l32)pu8Cur[l32Index] - (l32)pu8Nxt[l32Index];

        f32 f32SD1 = pf32Cur[l32Index] - pf32Pre[l32Index];
        f32 f32SD2 = pf32Cur[l32Index] - pf32Nxt[l32Index];

        pf32GD[l32Index] = /*sqrt*/((f32)ABS(l32GD1 * l32GD2));
        pf32SD[l32Index] = /*sqrt*/((f32)ABS(f32SD1 * f32SD2));
        f32GMean += pf32GD[l32Index];
        f32SMean += pf32SD[l32Index];       
    }
    f32GTh = 2.0f * f32GMean / (l32Height * l32Width);
    f32STh = 2.0f * f32SMean / (l32Height * l32Width);
    f32GTh = MAX(f32GTh, ptHandle->l32MinGrayDiff);
    f32STh = MAX(f32STh, ptHandle->l32MinStdDiff);

    for(l32Index = 0; l32Index < l32Height * l32Width; l32Index++)
    {
        if(pf32GD[l32Index] > f32GTh && pf32SD[l32Index] > f32STh) 
        {
            pu8Fg[l32Index] = FRAME_DIFF_FG;
        }
        else
        {
            pu8Fg[l32Index] = FRAME_DIFF_BG;
        }
    }
End:
    ptHandle->l32FrameNum++;
    return 0;
}
l32 FrameDiffClose(void *pvFrameDiffHandle)
{
    TFrameDiffHandle *ptHandle = NULL;
    if(NULL == pvFrameDiffHandle)
    {
        goto End;
    }
    ptHandle = (TFrameDiffHandle *)pvFrameDiffHandle;
    if(NULL != ptHandle->pu8ImgBuf)
    {
        free(ptHandle->pu8ImgBuf);
    }
    if(NULL != ptHandle->pf32StdBuf)
    {
        free(ptHandle->pf32StdBuf);
    }
    if(NULL != ptHandle->ppu8Img)
    {
        free(ptHandle->ppu8Img);
    }
    if(NULL != ptHandle->ppf32Std)
    {
        free(ptHandle->ppf32Std);
    }
    if(NULL != ptHandle->pvBuffer)
    {
        free(ptHandle->pvBuffer);
    }

    free(pvFrameDiffHandle);

End:
    return 0;
}

////////////================Integral image=====================/////////////
static void Integral(u8 *pu8Src, l32 *pl32Dst, l32 l32Width, l32 l32Hieght) 
{
    u8 *pu8Data = pu8Src;
    l32 l32SumRow;
    l32 l32Col, l32Row;

    // 1 row set zero
    pl32Dst += (l32Width + 1);

    // 2st row calculation
    for(l32Row = 1; l32Row < l32Hieght + 1; l32Row++, pl32Dst+=(l32Width + 1), pu8Data += l32Width)
    {
        l32SumRow = 0;
        for(l32Col = 1; l32Col < l32Width + 1; l32Col++)
        {
            l32SumRow += pu8Data[l32Col - 1];
            pl32Dst[l32Col] = pl32Dst[l32Col - l32Width - 1] + l32SumRow;
        }
    }
}

static void Integral2(u8 *pu8Src, l32 *pl32Dst, l32 l32Width, l32 l32Height) 
{
    u8 *pu8Data = pu8Src;
    l32 l32SumRow;
    l32 l32Col, l32Row;

    // 1st row set zero
    pl32Dst += (l32Width + 1);

    // 2nd row begins
    for(l32Row = 1; l32Row < l32Height + 1; l32Row++, pl32Dst += (l32Width + 1), pu8Data += l32Width)
    {
        l32SumRow = 0;
        for(l32Col = 1; l32Col < l32Width + 1; l32Col++)
        {
            l32SumRow += pu8Data[l32Col - 1] * pu8Data[l32Col - 1];
            pl32Dst[l32Col] = pl32Dst[l32Col - l32Width - 1] + l32SumRow;
        }
    }
}

//计算LC
static void CalcStd(u8 *pu8Gray, f32 *pf32Contrast, l32 * pl32IImg, l32 *pl32IImg2, l32 l32Width, l32 l32Height, l32 l32BlockW, l32 l32BlockH)
{
    l32 l32Row, l32Col, l32BlockSize;
    l32 l32Left, l32Right, l32Top, l32Bot;
    f32 f32Std, f32Var, f32EX, f32EX2;

    memset(pl32IImg, 0, (l32Width + 1) * (l32Height + 1) * sizeof(l32));
    memset(pl32IImg2, 0, (l32Width + 1) * (l32Height + 1) * sizeof(l32));

    //计算积分图像
    Integral(pu8Gray, pl32IImg, l32Width, l32Height);
    Integral2(pu8Gray, pl32IImg2, l32Width, l32Height);

    for(l32Row = 0; l32Row < l32Height; l32Row++)
    {
        for(l32Col = 0; l32Col < l32Width; l32Col++)
        {
            //compute local contrast
            l32Left  = MAX(l32Col - l32BlockW/2, 0);
            l32Right = MIN(l32Col + l32BlockW/2 + 1, l32Width);
            l32Top = MAX(l32Row - l32BlockH/2, 0);
            l32Bot = MIN(l32Row + l32BlockH/2 + 1, l32Height);
            l32BlockSize = (l32Right - l32Left) * (l32Bot - l32Top);

            f32EX = 1.0f * (pl32IImg[l32Bot * (l32Width + 1) + l32Right]  - pl32IImg[l32Bot * (l32Width + 1) + l32Left] 
            - pl32IImg[l32Top * (l32Width + 1) + l32Right] + pl32IImg[l32Top * (l32Width + 1) + l32Left]) / l32BlockSize;
            f32EX2 = 1.0f * (pl32IImg2[l32Bot * (l32Width + 1) + l32Right] - pl32IImg2[l32Bot * (l32Width + 1) + l32Left] 
            - pl32IImg2[l32Top * (l32Width + 1) + l32Right] + pl32IImg2[l32Top * (l32Width + 1) + l32Left]) / l32BlockSize;

            //计算方差
            f32Var = f32EX2 - f32EX * f32EX;
            f32Std = sqrt(f32Var);    
            *pf32Contrast = f32Std;
            pf32Contrast++;
        }
    }
}