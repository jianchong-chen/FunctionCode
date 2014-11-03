#ifndef _BG_DIFF_H_
#define _BG_DIFF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ai_defs.h"

#define DIFF_FG 255
#define DIFF_BG 0
#define BG_UPDATE_RATE 0.005
#define FG_UPDATE_RATE 0.0005

typedef struct
{
    l32 l32Height;
    l32 l32Width;
    l32 l32MinGrayDiff;
    l32 l32FrameNum;

    f32 *pf32Bg;
    f32 *pf32Buf;
}TBgDiffHandle;

l32 BgDiffOpen(void **ppvFrameDiffHandle, l32 l32Height, l32 l32Width, l32 l32MinGrayDiff = 5);
l32 BgDiffProcess(void *pvFrameDiffHandle, u8 *pu8Img, u8 *pu8Fg);
l32 BgDiffClose(void *pvFrameDiffHandle);

#ifdef __cplusplus
}
#endif
#endif