#ifndef _LBP_H__
#define _LBP_H__

#include "ai_defs.h"

l32 LBPModelOpen(void **pvvLBPHandle, l32 l32Width, l32 l32Height);
l32 LBPModelProcess(void *pvLBPHandle, u8 *pu8Img, u8 *pu8FG);
l32 LBPModelClose(void *pvLBPHandle);

#endif  //_LBP_H__

