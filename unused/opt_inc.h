#ifndef __OPT_INC_H_
#define __OPT_INC_H_

#include "ai_defs.h"
#include "uti_misc.h"

#ifdef __cplusplus

static void *operator new(size_t size)
{
	return AVMalloc(size);
}

static void * operator new[](size_t size)
{
	return AVMalloc(size);
}

static void operator delete(void *p)
{
	AVFree(p);
}

static void operator delete[](void *p)
{
	AVFree(p);
}
#endif

#ifdef ARCH_X86_WIN32
#include "intrinsic.h"

#define restrict 

//#define USE_OPENCV_DEBUG

extern "C" void UmpPrintf(const s8 *ps8Fmt, ...);
extern "C" u64 UmpGetHTime();

#endif

#ifdef USE_OPENCV_DEBUG
#include <opencv2/opencv.hpp>

#ifdef _DEBUG
#pragma comment(lib,"opencv_core247d.lib") 
#pragma comment(lib,"opencv_highgui247d.lib") 
#pragma comment(lib,"opencv_imgproc247d.lib") 
//#pragma comment(lib,"opencv_contrib247d.lib") 
//#pragma comment(lib,"opencv_calib3d247d.lib")
//#pragma comment(lib,"ipcsdk.lib")
//#pragma comment(lib,"virtual_pu.lib")
//#pragma comment(lib,"wsock32.lib")
#else
#pragma comment(lib,"opencv_core247.lib") 
#pragma comment(lib,"opencv_highgui247.lib") 
#pragma comment(lib,"opencv_imgproc247.lib") 
//#pragma comment(lib,"opencv_contrib247.lib") 
//#pragma comment(lib,"opencv_calib3d247.lib")
//#pragma comment(lib,"ipcsdk.lib")
//#pragma comment(lib,"virtual_pu.lib")
//#pragma comment(lib,"wsock32.lib")
#endif

static void ShowY8InOpenCV(const char * strFgName, u8 *pu8Src, l32 l32Width, l32 l32Height, l32 l32Stride, int waitKeyMs)
{
	int x, y;
	u8 * pu8Dst;
	IplImage *ImgY8 = cvCreateImage(cvSize(l32Width, l32Height), IPL_DEPTH_8U, 1);
	pu8Dst = (u8*)ImgY8->imageData;
	for(y=0; y<l32Height; y++)
	{
		memcpy(pu8Dst, pu8Src, l32Width);
		pu8Src += l32Stride;
		pu8Dst += ImgY8->widthStep;
	}
	cvShowImage(strFgName, ImgY8);
	if(waitKeyMs != 1) cvWaitKey(waitKeyMs);
	cvReleaseImage(&ImgY8);
}
#else
static void ShowY8InOpenCV(const char * strFgName, u8 *pu8Src, l32 l32Width, l32 l32Height, l32 l32Stride, int waitKeyMs)
{
}
#endif

#endif