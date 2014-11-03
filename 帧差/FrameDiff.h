#ifndef _FRAME_DIFF_H_
#define _FRAME_DIFF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ai_defs.h"

#define FRAME_DIFF_FG 255
#define FRAME_DIFF_BG 0

typedef struct
{
    l32 l32Height;
    l32 l32Width;
    l32 l32DiffFrameNum;
    l32 l32MinGrayDiff;
    l32 l32MinStdDiff;
    l32 l32FrameNum;

    u8 *pu8ImgBuf;
    f32 *pf32StdBuf;

    u8 **ppu8Img;     //»º´æµÄÍ¼Ïñ
    f32 **ppf32Std;   //»º´æµÄÍ¼Ïñ±ê×¼²î

    void *pvBuffer; // 2 * (H + 1) * (W + 1)
}TFrameDiffHandle;

//l32 FrameDiffOpen(void **ppvFrameDiffHandle, l32 l32Height, l32 l32Width);
l32 FrameDiffOpen(void **ppvFrameDiffHandle, l32 l32Height, l32 l32Width, l32 l32DiffFrameNum = 4, l32 l32MinGrayDiff = 5, l32 l32MinStdDiff = 2);
l32 FrameDiffProcess(void *pvFrameDiffHandle, u8 *pu8Img, u8 *pu8Fg);
l32 FrameDiffClose(void *pvFrameDiffHandle);

#ifdef __cplusplus
}
#endif
#endif