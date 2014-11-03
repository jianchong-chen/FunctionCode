#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "opt_inc.h"
#include "lcgmm.h"


#if defined (ARCH_C667X)
#define _UNPKBU4(u32DataC) _unpkbu4(u32DataC)
#define _DADD2(u64Src1, u64Src2) _dadd2(u64Src1, u64Src2)
#define _DSUB2(u64Src1, u64Src2) _dsub2(u64Src1, u64Src2)
#define _DSHRU(u64Src1, u32Src2) _dshru(u64Src1, u32Src2)
#define _DSHL(u64Src1, u32Src2) _dshl(u64Src1, u32Src2)
#define _DADD(u64Src1, u64Src2) _dadd(u64Src1, u64Src2)
#define _DPACKHL2(u64Src1,u64Src2) _dpackhl2(u64Src1,u64Src2)
#define _UNPKHU2(u32Src1) _unpkhu2(u32Src1)
#define _MPYU2(u32Src1, u32Src2) _mpyu2(u32Src1, u32Src2)
#define _DCMPGTU4(s64src1, s64src2) _dcmpgtu4(s64src1, s64src2)
#define _DXPND4(u32src) _dxpnd4(u32src)
#endif

#if defined (ARCH_C674X)
#define _UNPKBU4(u32DataC) _itoll(_unpkhu4(u32DataC), _unpklu4(u32DataC))
#define _DADD2(u64Src1, u64Src2) _itoll(_add2(_hill(u64Src1), _hill(u64Src2)), _add2(_loll(u64Src1), _loll(u64Src2)))
#define _DSUB2(u64Src1, u64Src2) _itoll(_sub2(_hill(u64Src1), _hill(u64Src2)), _sub2(_loll(u64Src1), _loll(u64Src2)))
#define _DSHRU(u64Src1, u32Src2) _itoll((_hill(u64Src1) >> u32Src2), (_loll(u64Src1) >> u32Src2))
#define _DSHL(u64Src1, u32Src2) _itoll((_hill(u64Src1) << u32Src2), (_loll(u64Src1) << u32Src2))
#define _DADD(u64Src1, u64Src2) _itoll(_hill(u64Src1) + _hill(u64Src2), _loll(u64Src1) + _loll(u64Src2))
#define _DPACKHL2(u64Src1,u64Src2) _itoll(_packhl2(_hill(u64Src1), _hill(u64Src2)), _packhl2(_loll(u64Src1), _loll(u64Src2)))
#define _UNPKHU2(u32Src1) _itoll((u32Src1>>16), (u32Src1 & 0xFFFF))
#define _MPYU2(u32Src1, u32Src2) _itoll(_mpyhu(u32Src1, u32Src2),_mpyu(u32Src1, u32Src2))
#define _DCMPGTU4(s64src1, s64src2) (_cmpgtu4(_hill(s64src1), _hill(s64src2)) << 4)|(_cmpgtu4(_loll(s64src1), _loll(s64src2)))
#define _DXPND4(u32src) _itoll(_xpnd4(u32src >> 4), _xpnd4(u32src))
#endif


void iimg_4_n(u8 * restrict pu8Src, s16 *ps16Dst, u32 *pu32Dst, int imageHeight, int imageWidth) 
{
    u64 u64DataC, u64DataP1;
    u64 u64DataP2, u64DataP3;

    u64 u64Res1, u64Res2;
    u64 u64Res3, u64Res4;

    u64 u64Data1, u64Data2;
    u64 u64Data3, u64Data4;

    u32 u32DataC, u32DataP1;
    u32 u32DataP2, u32DataP3;
    u32 u32Left, u32Right; 

    l32 l32I, l32J;
    l32 l32SrcWidth, l32DstWidth;
    u8 *pu8SrcC, *pu8SrcP1;
    u8 *pu8SrcP2, *pu8SrcP3;
    s16 *ps16DstC;
    u32 *pu32DstC;

    s16 s16Sum;

    pu8SrcC = pu8Src;
    pu8SrcP1 = pu8Src + imageWidth;
    pu8SrcP2 = pu8SrcP1 + imageWidth;
    pu8SrcP3 = pu8SrcP2 + imageWidth;

    ps16DstC = ps16Dst;
    pu32DstC = pu32Dst;

    l32SrcWidth = imageWidth << 2;
    l32DstWidth = imageWidth >> 2;

    for (l32I = 0; l32I <= imageHeight - 4; l32I += 4)
    {
#if defined (ARCH_C667X) || defined (ARCH_C674X)
        for (l32J = 0; l32J <= imageWidth - 4; l32J += 4)
        {
            u32DataC = _mem4(pu8SrcC + l32J);
            u64DataC = _UNPKBU4(u32DataC);

            u32DataP1 = _mem4(pu8SrcP1 + l32J);
            u64DataP1 = _UNPKBU4(u32DataP1);

            u32DataP2 = _mem4(pu8SrcP2 + l32J);
            u64DataP2 = _UNPKBU4(u32DataP2);

            u32DataP3 = _mem4(pu8SrcP3 + l32J);
            u64DataP3 = _UNPKBU4(u32DataP3);

            u64Res1 = _DADD2(u64DataC, u64DataP1);
            u64Res2 = _DADD2(u64DataP2, u64DataP3);
            u64Res3 = _DADD2(u64Res1, u64Res2);

            u32Left = _loll(u64Res3);
            u32Right = _hill(u64Res3);

            u32Left = _add2(u32Left, u32Right);

            s16Sum = _ext(u32Left, 16, 16) + _ext(u32Left, 0, 16);

            ps16DstC[l32J >> 2] = s16Sum;

            u64Res1 = _mpyu4ll(u32DataC, u32DataC);
            u64Res2 = _mpyu4ll(u32DataP1, u32DataP1);
            u64Res3 = _mpyu4ll(u32DataP2, u32DataP2);
            u64Res4 = _mpyu4ll(u32DataP3, u32DataP3);

            u64Data1 = _UNPKHU2(_loll(u64Res1));
            u64DataC = _UNPKHU2(_hill(u64Res1));

            u64Data2 = _UNPKHU2(_loll(u64Res2));
            u64DataP1 = _UNPKHU2(_hill(u64Res2));

            u64Data3 = _UNPKHU2(_loll(u64Res3));
            u64DataP2 = _UNPKHU2(_hill(u64Res3));

            u64Data4 = _UNPKHU2(_loll(u64Res4));
            u64DataP3 = _UNPKHU2(_hill(u64Res4));

            u64Data1 = _DADD(u64Data1, u64Data2);
            u64Data2 = _DADD(u64Data3, u64Data4);
            u64Data1 = _DADD(u64Data1, u64Data2);

            u64DataC = _DADD(u64DataC, u64DataP1);
            u64DataP1 = _DADD(u64DataP2, u64DataP3);
            u64DataC = _DADD(u64DataC, u64DataP1);

            u64DataC = _DADD(u64Data1, u64DataC);

            u32Left = _loll(u64DataC);
            u32Right = _hill(u64DataC);

            u32Left = u32Left + u32Right;

            pu32DstC[l32J >> 2] = u32Left;
        }
#else
        l32 i2;
        u8 *pu8Temp;

        memset(ps16DstC, 0, sizeof(s16)*imageWidth/4);
        memset(pu32DstC, 0, sizeof(u32)*imageWidth/4);


        pu8Temp = pu8SrcC;
        for(i2 = 0; i2 < 4; i2++)
        {
            for (l32J = 0; l32J <= imageWidth - 4; l32J += 4)
            {
                s16Sum = pu8Temp[l32J] + pu8Temp[l32J+1] + pu8Temp[l32J+2] + pu8Temp[l32J+3];
                ps16DstC[l32J >> 2] += s16Sum;

                u32Left = pu8Temp[l32J]*pu8Temp[l32J] 
                        + pu8Temp[l32J+1]*pu8Temp[l32J+1]
                        + pu8Temp[l32J+2]*pu8Temp[l32J+2]
                        + pu8Temp[l32J+3]*pu8Temp[l32J+3];

                pu32DstC[l32J >> 2] += u32Left;
            }
            pu8Temp += imageWidth;
        }

#endif

        pu8SrcC += l32SrcWidth;
        pu8SrcP1 += l32SrcWidth;
        pu8SrcP2 += l32SrcWidth;
        pu8SrcP3 += l32SrcWidth;
        ps16DstC += l32DstWidth;
        pu32DstC += l32DstWidth;
    }
}

void update_bg_img(u8 *pu8YUV, u8 *pu8MeanBg, u8 *pu8Mask, u16 * restrict pu16MeanBg, 
                   int l32Width, int l32Height,
                   float f32LcLearnRate)
{
    u8 * restrict pu8MeanBgLc = pu8MeanBg;
    s8 *  restrict pOutBuffer = (s8*)pu8Mask; // pContext->pOutBuffer;  //用s8*可以容易扩展符号位
    u8 * restrict pu8GmmPlusGray = pu8YUV; // pContext->pYuvBuffer;

    int i, iMask, iMask0, iMask1;
    u32 u32Fg0, u32Fg1;
    u32 u32MeanBg0, u32MeanBg1;
    u32 u32Mask, u32Mask0, u32Mask1;
    u32 u32ConstW1, u32ConstW2;
    u64 u64Data0, u64Data1, u64Data2, u64Data3;
    u64 u64MeanBg, u64Bg, u64Mask;
    u32 u32MeanBg, u32Bg;
    u32 u32LearnRate = (u32)(f32LcLearnRate * 128);

	static int bFirstFrame = 1;
	if(bFirstFrame)
	{
		bFirstFrame = 0;
		for(i = 0; i < (l32Width * l32Height); i++)
		{
			pu16MeanBg[i] = pu8GmmPlusGray[i] << 8;
		}
	}
#if defined (ARCH_C667X)
    /*
        假设原始均值背景为M,定点化倍率为S,现在的输入图像发生变化x(或正或负)，则M产生更新如下：
            M = (M*w1 + (x*S)*w2)/(w1+w2)
              = M + (x*S-M)*w2/(w1+w2);
        因为u8图像，x是整数，因此为了让M能够产生+1的变化，有 (x*S) > M + ((w1+w2)/w2)
        换言之就是在IIR背景更新过程中，如果外部x不变，则M会逐步更新，但是只能更新到距离x尚有一点差别就会停止了
        这个差别就是((w1+w2)/w2)/S，S取值越大，这个差别就越小。
                例如归一化权值w2=0.1的话，S取1，也就是完全用整数计算的话，这个差别达到9，因为小于等于9的差别
                无法产生进位效果。比方说M为0，输入x为9，上面的更新过程总是不能使M增加
                
                但是S取10时，这个差别就是1，仍以M为0，输入x为1为例：
                    第一轮更新   M = 0 + (1*10-0)*0.1; = 1
                    第二轮更新   M = 1 + (1*10-1)*0.1; = 1   无法进一步更新

                也就是一旦(x*S-M)<(w1+w2)/w2就停止更新

        要更新到M的整数部分完全跟x*S一致，就需要M跟x*S的误差在正负(S/2)以内(这样四舍五入就可以达到一致目的)
        也就是(w1+w2)/(w2)这个停止更新的阈值要小于S/2。也就是S > (2*(w1+w2)/w2)
    */
    u32ConstW1 = _pack2(128 - u32LearnRate,128 - u32LearnRate);
    u32ConstW2 = _pack2(u32LearnRate*256, u32LearnRate*256);
    for(i = 0; i < (l32Width * l32Height); i+=4)
    {
        //u32MeanBg的精度是定点数,15bit小数位
        
        //读入四个背景点，四个输入点，16bit无符号数
        u64Bg = _amem8(pu16MeanBg + i);
        u32Fg1 = _amem4(pu8GmmPlusGray + i);
        u32Fg0 = _unpklu4(u32Fg1);
        u32Fg1 = _unpkhu4(u32Fg1);
        
        //加权滤波
        u64Data0 = _mpyu2(_loll(u64Bg), u32ConstW1);
        u64Data1 = _mpyu2(_hill(u64Bg), u32ConstW1);
        u64Data2 = _mpyu2(u32Fg0, u32ConstW2);
        u64Data3 = _mpyu2(u32Fg1, u32ConstW2);
        u64Data0 = _dadd(u64Data0, u64Data2);
        u64Data1 = _dadd(u64Data1, u64Data3);

        u64Data0 = _dshru(u64Data0, 7);
        u64Data1 = _dshru(u64Data1, 7);

        
        u64MeanBg = _itoll(_pack2(_hill(u64Data1), _loll(u64Data1)),
                            _pack2(_hill(u64Data0), _loll(u64Data0)));

        u32Mask = _amem4(pOutBuffer + i);
        u64Mask = _unpkbu4(u32Mask);
        u64Mask = _dshl(u64Mask, 8);
        u64Mask = _dshr2(u64Mask, 8);

        //合并得到最终结果
        u64MeanBg = (u64MeanBg & (~u64Mask)) | (u64Bg & u64Mask);

        //保存回去
        _amem8(pu16MeanBg + i) = u64MeanBg;
        _amem4(pu8MeanBgLc + i) = _packh4(_hill(u64MeanBg), _loll(u64MeanBg));
        
        /*

        u32Bg = pu16MeanBg[i];
        u32MeanBg = (u32Bg * (256-LC_LEARN_RATE) + pu8GmmPlusGray[i] * (LC_LEARN_RATE*(1<<8))) >> 8;
        iMask = pOutBuffer[i];

        u32MeanBg = (u32Bg & iMask) + (u32MeanBg & (~iMask));

        pu16MeanBg[i] = u32MeanBg;
        pu8MeanBg[i] = u32MeanBg >> 8;
        */
    }
#elif defined (ARCH_C674X)
    u32ConstW1 = _pack2(128 - u32LearnRate,128 - u32LearnRate);
    u32ConstW2 = _pack2(u32LearnRate*256, u32LearnRate*256);
    for(i = 0; i < (l32Width * l32Height); i+=2)
    {
        //读入两个背景点，两个输入点，16bit无符号数
        u32Bg = _amem4(pu16MeanBg + i);
        u32Fg1 = _amem2(pu8GmmPlusGray + i);
        u32Fg0 = _unpklu4(u32Fg1);
        
        //加权滤波
        u64Data0 = _MPYU2(u32Bg, u32ConstW1);
        u64Data2 = _MPYU2(u32Fg0, u32ConstW2);
        u64Data0 = _DADD(u64Data0, u64Data2);
        u64Data0 = _DSHRU(u64Data0, 7);

        u32MeanBg = _pack2(_hill(u64Data0), _loll(u64Data0));

        u32Mask = _amem2(pOutBuffer + i);
        u32Mask = _unpklu4(u32Mask);
        u32Mask = u32Mask << 8;
        u32Mask = _shr2(u32Mask, 8);

        //合并得到最终结果
        u32MeanBg = (u32MeanBg & (~u32Mask)) | (u32Bg & u32Mask);

        //保存回去
        _amem4(pu16MeanBg + i) = u32MeanBg;
        _amem2(pu8MeanBgLc + i) = _packh4(0, _loll(u32MeanBg));
    }
#else
    for(i = 0; i < (l32Width * l32Height); i++)
    {
        if(pu8Mask[i] == 0)
        {
            u32MeanBg0 = (pu16MeanBg[i] * (128-u32LearnRate) + pu8GmmPlusGray[i] * (u32LearnRate*(256))) >> 7;
            pu16MeanBg[i] = u32MeanBg0;
        }
        else
        {
            u32MeanBg0 = pu16MeanBg[i];
        }
        pu8MeanBgLc[i] = u32MeanBg0 >> 8;
    }
#endif
}

static INLINE void Sobel1Line( u8 * restrict pu8Src0,
                               u8 * restrict pu8Src1,
                               u8 * restrict pu8Src2,
                               l32 l32Width,
                               l32 * pl32SobelHV)
{
    l32 l32x, l32x2, l32tmp0, l32tmp1, l32tmp2, l32H, l32V;
    u64 u64DataPre_Left, u64DataPre_Right, u64DataCur_Left, u64DataCur_Right,u64DataNex_Left, u64DataNex_Right;
    u64 u64Pre_Left0, u64Pre_Left1, u64Cur_Left0, u64Cur_Left1, u64Nex_Left0, u64Nex_Left1;
    u64 u64Pre_Right0, u64Pre_Right1, u64Cur_Right0, u64Cur_Right1, u64Nex_Right0, u64Nex_Right1;
    u64 u64Pre_V0, u64Pre_V1, u64Cur_V0, u64Cur_V1, u64Nex_V0, u64Nex_V1;
    u64 u64Temp0, u64Temp1, u64Temp2, u64Temp3;
    u64 u64Top, u64Bot, u64Top0, u64Top1, u64Bot0, u64Bot1;
    u64 u64V0, u64V1, u64V2, u64V3;
    u64 u64H0, u64H1, u64H2, u64H3, u64H4, u64H5;

    l32x = 1;
#if (defined (ARCH_C667X) || defined (ARCH_C674X))
    for(; l32x < l32Width - 8; l32x+=8)
    {
		u64DataPre_Left = _mem8_const(pu8Src0 + l32x - 1); //a7 a6 a5 a4 a3 a2 a1 a0
		u64DataPre_Right = _mem8_const(pu8Src0 + l32x + 1);//a9 a8 a7 a6 a5 a4 a3 a2
		u64DataCur_Left = _mem8_const(pu8Src1 + l32x - 1); //b7 b6 b5 b4 b3 b2 b1 b0
		u64DataCur_Right = _mem8_const(pu8Src1 + l32x + 1);//b9 b8 b7 b6 b5 b4 b3 b2
		u64DataNex_Left = _mem8_const(pu8Src2 + l32x - 1); //c7 c6 c5 c4 c3 c2 c1 c0
		u64DataNex_Right = _mem8_const(pu8Src2 + l32x + 1);//c9 c8 c7 c6 c5 c4 c3 c2

		u64Pre_Left0 = _UNPKBU4(_loll(u64DataPre_Left));//0 a3 0 a2 0 a1 0 a0
		u64Pre_Left1 = _UNPKBU4(_hill(u64DataPre_Left));//0 a7 0 a6 0 a5 0 a4

		u64Pre_Right0 = _UNPKBU4(_loll(u64DataPre_Right));//0 a5 0 a4 0 a3 0 a2
		u64Pre_Right1 = _UNPKBU4(_hill(u64DataPre_Right));//0 a9 0 a8 0 a7 0 a6

		u64Cur_Left0 = _UNPKBU4(_loll(u64DataCur_Left));//0 b3 0 b2 0 b1 0 b0
		u64Cur_Left1 = _UNPKBU4(_hill(u64DataCur_Left));//0 b7 0 b6 0 b5 0 b4

		u64Cur_Right0 =  _UNPKBU4(_loll(u64DataCur_Right));//0 b5 0 b4 0 b3 0 b2
		u64Cur_Right1 =  _UNPKBU4(_hill(u64DataCur_Right));//0 b9 0 b8 0 b7 0 b6

		u64Nex_Left0 = _UNPKBU4(_loll(u64DataNex_Left));//0 c3 0 c2 0 c1 0 c0
		u64Nex_Left1 = _UNPKBU4(_hill(u64DataNex_Left));//0 c7 0 c6 0 c5 0 c4

		u64Nex_Right0 =  _UNPKBU4(_loll(u64DataNex_Right));//0 c5 0 c4 0 c3 0 c2
		u64Nex_Right1 =  _UNPKBU4(_hill(u64DataNex_Right));//0 c9 0 c8 0 c7 0 c6


		u64Pre_V0 = _DSUB2(u64Pre_Right0, u64Pre_Left0);//a2 - a0  a3 - a1
		u64Pre_V1 = _DSUB2(u64Pre_Right1, u64Pre_Left1);

		u64Cur_V0 = _DSUB2(u64Cur_Right0, u64Cur_Left0);//b2 - b0  b3 - b1
		u64Cur_V1 = _DSUB2(u64Cur_Right1, u64Cur_Left1);

		u64Nex_V0 = _DSUB2(u64Nex_Right0, u64Nex_Left0);//c2 - c0  c3 - c1
		u64Nex_V1 = _DSUB2(u64Nex_Right1, u64Nex_Left1);

		u64V0 = _DADD2(u64Pre_V0, u64Nex_V0);
		u64V1 = _DADD2(u64Pre_V1, u64Nex_V1);

		u64V0 = _DADD2(u64V0, u64Cur_V0);
		u64V1 = _DADD2(u64V1, u64Cur_V1);

		u64V0 = _DADD2(u64V0, u64Cur_V0);
		u64V1 = _DADD2(u64V1, u64Cur_V1);

		u64Temp0 = _DSHRU(u64DataPre_Left, 8); //0 a7 a6 a5 0 a3 a2 a1
		u64Temp1 = _DSHL(u64DataPre_Right, 8);//a8 a7 a6 0 a4 a3 a2 0

		u64Temp2 = _DSHRU(u64DataNex_Left, 8);
		u64Temp3 = _DSHL(u64DataNex_Right, 8);

		u64Top = _DPACKHL2(u64Temp1, u64Temp0);//a8 a7 a6 a5 a4 a3 a2 a1
		u64Bot = _DPACKHL2(u64Temp3, u64Temp2);//c8 c7 c6 c5 c4 c3 c2 c1

		u64Top0 = _UNPKBU4(_loll(u64Top));
		u64Top1 = _UNPKBU4(_hill(u64Top));

		u64Bot0 = _UNPKBU4(_loll(u64Bot));
		u64Bot1 = _UNPKBU4(_hill(u64Bot));

		u64H0 = _DSUB2(u64Nex_Left0, u64Pre_Left0);  //c0 - a0  c1 - a1  c2 - a2  c3 - a3
		u64H1 = _DSUB2(u64Nex_Left1, u64Pre_Left1);  //c4 - a4  c5 - a5  c6 - a6  c7 - a7
		u64H2 = _DSUB2(u64Bot0, u64Top0);            //c1 - a1  c2 - a2  c3 - a3  c4 - a4
		u64H3 = _DSUB2(u64Bot1, u64Top1);            //c5 - a5  c6 - a6  c7 - a7  c8 - a8
		u64H4 = _DSUB2(u64Nex_Right0, u64Pre_Right0);//c2 - a2  c3 - a3  c4 - a4  c5 - a5
		u64H5 = _DSUB2(u64Nex_Right1, u64Pre_Right1);//c6 - a6  c7 - a7  c8 - a8  c9 - a9

		u64H0 = _DADD2(u64H0, u64H2);
		u64H1 = _DADD2(u64H1, u64H3);
		u64H0 = _DADD2(u64H0, u64H2);
		u64H1 = _DADD2(u64H1, u64H3);
		u64H0 = _DADD2(u64H0, u64H4);
		u64H1 = _DADD2(u64H1, u64H5);

        //HV交替存放
        _amem8((void *)(pl32SobelHV + l32x)) = _dpack2(_loll(u64H0), _loll(u64V0));
        _amem8((void *)(pl32SobelHV + l32x + 2)) = _dpack2(_hill(u64H0), _hill(u64V0));

        _amem8((void *)(pl32SobelHV + l32x + 4)) = _dpack2(_loll(u64H1), _loll(u64V1));
        _amem8((void *)(pl32SobelHV + l32x + 6)) = _dpack2(_hill(u64H1), _hill(u64V1));

    }

#endif
    //求一行的sobel H和V数据
    for(; l32x < l32Width - 1; l32x++)
    {
        l32tmp0 = (pu8Src0[l32x-1] - pu8Src0[l32x+1]);
        l32tmp1 = (pu8Src1[l32x-1] - pu8Src1[l32x+1]);
        l32tmp2 = (pu8Src2[l32x-1] - pu8Src2[l32x+1]);
        l32V = (l32tmp0 + l32tmp1*2 + l32tmp2);

        l32tmp0 = (pu8Src0[l32x-1] - pu8Src2[l32x-1]);
        l32tmp1 = (pu8Src0[l32x] - pu8Src2[l32x]);
        l32tmp2 = (pu8Src0[l32x+1] - pu8Src2[l32x+1]);
        l32H = (l32tmp0 + l32tmp1*2 + l32tmp2);

        pl32SobelHV[l32x] = (l32V & 0xFFFF) + (l32H << 16);

        //ps16SobelV[l32x] = l32i0;
        //ps16SobelH[l32x] = l32i2;
    }
}

static INLINE void DirCompare1Line( l32 * restrict pl32SrcHV,
								 l32 * restrict pl32BgHV,
								 u8 * restrict pu8Dst, l32 l32Width,
								 int TgFix, int TsFix)
{
	l32 l32x;
    int HV0, HV1, dotp, mod_cur, mod_bg, not_flat, vote;
    int dotp1, mod_cur1, mod_bg1, not_flat1, vote1;
    u64 u64HV0, u64HV1, u64Temp0;
    u32 TgFix2 = _pack2(TgFix, TgFix);
    u32 TsFix2 = _pack2(TsFix, TsFix);

//#if defined (ARCH_C667X)
	for(l32x = 0; l32x < l32Width; l32x++)
	{
		//对比本行的对应点的sobel
		/*
		int V = pl32SrcHV[l32x] >> 16;
		int H = (pl32SrcHV[l32x]<<16)>>16;
		int Vbg = pl32BgHV[l32x] >> 16;
		int Hbg = (pl32BgHV[l32x]<<16)>>16;
		dotp    = V*Vbg  + H*Hbg;
		mod_cur = V*V    + H*H;
		mod_bg = Vbg*Vbg + Hbg*Hbg;
		*/
		HV0 = pl32SrcHV[l32x];
		HV1 = pl32BgHV[l32x];
		dotp    = _dotp2(HV0, HV1);
		mod_cur = _dotp2(HV0, HV0);
		mod_bg  = _dotp2(HV1, HV1);

		not_flat = (mod_cur) > TgFix;      //非平坦点(平坦点取背景建模给出的结果)
		//int vote = ((2*dotp) <= ((mod_cur + mod_bg)>>6) * TsFix);//点积不够大，说明主方向不一致，投1票给前景
		vote = (dotp < ((mod_cur + mod_bg)>>6) * TsFix);//点积不够大，说明主方向不一致，投1票给前景
		/*
			方向矢量v,u
			a=|v|;  a是v的幅度
			b=|u|;  b是u的幅度
			q = v,u夹角
			2*v*u/(|v|^2 + |u|^2)代表什么呢？
			2*(v.u)/(|v|^2 + |u|^2) = 2abcos(q)/(a^2 + b^2)

			此值在q=0时取得最大值，而且在a=b时取得最大。
			q=0说明u,v方向一致，a=b说明sobel边缘强度没有变化

			因此不论如何，此值越大说明越可能是背景
		*/
		//vote = ((float)dotp*(float)dotp) <= (((float)(mod_cur)*float(mod_bg))*0.64f);//点积不够大，说明主方向不一致，投1票给前景
		pu8Dst[l32x] |= vote & not_flat;
	}
}
static INLINE void DensityCompare1Line(u8 * restrict pu8Src1, 
									   u8 * restrict pu8Bg1, 
									   u8 * restrict pu8Dst, l32 l32Width, int TdFix)
{
	int l32x;
	u32 u32Data0, u32Data1, u32Data2, u32Data3;
	u32 u32TdFix;

	u32TdFix = TdFix | (TdFix << 16);
	u32TdFix = _packl4(u32TdFix, u32TdFix);

	//当前图像与背景图像差值的绝对值是否大于阈值TdFix
	for(l32x = 0; l32x < l32Width; l32x+= 4)
	{
		u32Data0 = *(u32 *)(pu8Src1 + l32x);
		u32Data1 = *(u32 *)(pu8Bg1 + l32x);
		u32Data2 = _subabs4(u32Data0, u32Data1);
		u32Data3 = _cmpgtu4(u32Data2, u32TdFix);
		u32Data3 = _xpnd4(u32Data3);
		*(u32 *)(pu8Dst + l32x) |= (u32Data3 & 0x01010101);
	}
}

//对输入两幅图像的sobel相似度比较
void SobelCompare_n(u8 *pu8Src, u8 *pu8Bg, l32 l32Width, l32 l32Height, l32 l32Stride, 
				  int TgFix, int TsFix, int TdFix, u8 *pu8Dst, s16 * ps16Temp4Line, BOOL bClear)
{
	l32 l32x, l32y;

	u8 * restrict pu8Src0;
	u8 * restrict pu8Src1;
	u8 * restrict pu8Src2;

	u8 * restrict pu8Bg0;
	u8 * restrict pu8Bg1;
	u8 * restrict pu8Bg2;

	l32 * restrict pl32SrcHV = (l32*)ps16Temp4Line;
	l32 * restrict pl32BgHV = pl32SrcHV + l32Width;

	pu8Src1 = pu8Src;
	pu8Bg1 = pu8Bg;
	for(l32y = 0; l32y < l32Height; 
		l32y++, 
		pu8Src1 += l32Width,
		pu8Bg1 += l32Width,
		pu8Dst += l32Width)
	{
		pu8Src0 = pu8Src1;
		pu8Bg0 = pu8Bg1;
		if(l32y > 0) 
		{
			pu8Src0 -= l32Width;
			pu8Bg0 -= l32Width;
		}

		pu8Src2 = pu8Src1;
		pu8Bg2 = pu8Bg1;
		if(l32y < l32Height-1)
		{
			pu8Src2 += l32Width;
			pu8Bg2 += l32Width;
		}

		//求一行的sobel H和V数据
		Sobel1Line(pu8Src0, pu8Src1, pu8Src2, l32Width, pl32SrcHV);
		Sobel1Line(pu8Bg0, pu8Bg1, pu8Bg2, l32Width,    pl32BgHV);
		if(bClear)
			memset(pu8Dst, 0, l32Width);

		DirCompare1Line(pl32SrcHV, pl32BgHV, pu8Dst, l32Width, TgFix,TsFix);

		DensityCompare1Line(pu8Src1, pu8Bg1, pu8Dst, l32Width, TdFix);
	}
}

static void Sum5x5(u8 *pu8Src, l32 l32Width, l32 l32Height, u8 *pu8Dst, u8 *pu8Temp5Line)
{
    l32 l32x, l32y, l32id, l32Sum;
    u8 * restrict pu8HResult;
    u8 * restrict pu8Line0;
    u8 * restrict pu8Line1;
    u8 * restrict pu8Line2;
    u8 * restrict pu8Line3;
    u8 * restrict pu8Line4;
    u64 u64Data0, u64Data1, u64Data2, u64Data3, u64Data4;
    u32 u32DataL0, u32DataH0;
    u32 u32Data0, u32Data1, u32Data2, u32Data3, u32Data4;

    memset(pu8Temp5Line, 0, l32Width*5);
    pu8Line0 = pu8Temp5Line;
    pu8Line1 = pu8Line0 + l32Width;
    pu8Line2 = pu8Line1 + l32Width;
    pu8Line3 = pu8Line2 + l32Width;
    pu8Line4 = pu8Line3 + l32Width;

    for(l32y = 0, l32id = 0; l32y < l32Height + 2; l32y++, pu8Src+= l32Width)
    {
        pu8HResult = pu8Temp5Line + (l32id) * l32Width;
        l32id ++;
        if(l32id > 4) l32id = 0;

        if(l32y < l32Height)
        {
            pu8HResult[0] = pu8HResult[1] = 0;
            l32x = 2;
#if (defined (ARCH_C667X) || defined (ARCH_C674X))
            //当前行水平5点求和
            for(; (l32x + 8) < l32Width; l32x += 8)
            {
                u64Data0 = _mem8_const(pu8Src + l32x - 2);
                u64Data1 = _mem8_const(pu8Src + l32x - 1);
                u64Data2 = _mem8_const(pu8Src + l32x);
                u64Data3 = _mem8_const(pu8Src + l32x + 1);
                u64Data4 = _mem8_const(pu8Src + l32x + 2);

                u32DataL0 = _add4(_add4(_loll(u64Data0), _loll(u64Data1)), _add4(_loll(u64Data2), _loll(u64Data3)));
                u32DataL0 = _add4(u32DataL0, _loll(u64Data4));

                u32DataH0 = _add4(_add4(_hill(u64Data0), _hill(u64Data1)), _add4(_hill(u64Data2), _hill(u64Data3)));
                u32DataH0 = _add4(u32DataH0, _hill(u64Data4));

                _mem8(pu8HResult + l32x) = _itoll(u32DataH0, u32DataL0);
            }
#endif
            if(l32x < l32Width)
            {
                l32Sum = pu8Src[l32x - 2] + pu8Src[l32x - 1] + pu8Src[l32x] + pu8Src[l32x + 1] + pu8Src[l32x + 2];
                for(; l32x < l32Width; l32x ++)
                {
                    pu8HResult[l32x] = l32Sum;
                    l32Sum -= pu8Src[l32x - 2];
                    l32Sum += pu8Src[l32x + 3];
                }
            }
        }
        else
        {
            //超出画面部分的sobel比较结果使用0代替（0代表跟背景一致）
            memset(pu8HResult, 0, l32Width);
        }

        //垂直5行求和
        if(l32y > 1)
        {
            //数据积累满5行（需要超前2行） -2,-1, 0, 1, 2
            l32x = 0;
#if (defined (ARCH_C667X) || defined (ARCH_C674X))
            for(; l32x + 8< l32Width; l32x += 8)
            {
                u64Data0 = _amem8_const(pu8Line0 + l32x);
                u64Data1 = _amem8_const(pu8Line1 + l32x);
                u64Data2 = _amem8_const(pu8Line2 + l32x);
                u64Data3 = _amem8_const(pu8Line3 + l32x);
                u64Data4 = _amem8_const(pu8Line4 + l32x);

                u32DataL0 = _add4(_add4(_loll(u64Data0), _loll(u64Data1)), _add4(_loll(u64Data2), _loll(u64Data3)));
                u32DataL0 = _add4(u32DataL0, _loll(u64Data4));

                u32DataH0 = _add4(_add4(_hill(u64Data0), _hill(u64Data1)), _add4(_hill(u64Data2), _hill(u64Data3)));
                u32DataH0 = _add4(u32DataH0, _hill(u64Data4));

                _amem8(pu8Dst + l32x) = _itoll(u32DataH0, u32DataL0);
            }
#endif
            if(l32x < l32Width)
            {
                for(; l32x < l32Width; l32x ++)
                {
                    pu8Dst[l32x] = pu8Line0[l32x]
                    + pu8Line1[l32x]
                    + pu8Line2[l32x]
                    + pu8Line3[l32x]
                    + pu8Line4[l32x];
                }
            }
            pu8Dst += l32Width;
        }
    }
}

static void ClearFalseFg(u8 * restrict pu8Fg, u8 * restrict pu8VoteMap, l32 l32Th, l32 l32Width, l32 l32Height)
{
    l32 l32x, l32y;
    u64 u64Data0, u64Data1, u64Data2, u64Th;
    u32 u32Data0, u32Data1;

    u32Data0 = (l32Th << 8) | l32Th;
    u32Data0 = (u32Data0 << 16) | u32Data0;
    u64Th = _itoll(u32Data0, u32Data0);

    for(l32y = 0; l32y < l32Height; l32y++)
    {
        l32x = 0;
#if (defined (ARCH_C667X) || defined (ARCH_C674X))
        for(; l32x < l32Width-7; l32x += 8)
        {
            u64Data0 = _amem8_const(pu8Fg + l32x);
            u64Data1 = _amem8_const(pu8VoteMap + l32x);

            //大于此阈值的得以保留原有前景结果
            u32Data0 = _cmpgtu4(_loll(u64Data1), _loll(u64Th));
            u32Data1 = _cmpgtu4(_hill(u64Data1), _hill(u64Th));
            u64Data2 = _itoll(_xpnd4(u32Data1), _xpnd4(u32Data0));

            //_amem8(pu8Fg + l32x) = (u64Data0 & u64Data2);
            _amem8(pu8Fg + l32x) = u64Data2;
        }
#endif
		//有疑问，会造成ＧＭＭ＋前景甚至多于ＧＭＭ结果
        for(; l32x < l32Width; l32x ++)
        {
            if((pu8VoteMap[l32x] <= l32Th))
            {
                pu8Fg[l32x] = 0;
            }
            else
            {
                pu8Fg[l32x] = 255;
            }			
        }
        pu8Fg += l32Width;
        pu8VoteMap += l32Width;
    }
}

void CLcgmm::gmmplus(u8 *pu8YUV, u8 *pu8MeanBg, u8 *pu8Mask)
{
    l32 i, j, m;
    l32 l32I, l32J, l32K;
    l32 nBlobPix, l32IimgWidth;
    u64 u64Data;
    u8 * restrict pu8OutBuffer;
    u8 * restrict pu8OutBufferBase;
    u8 * restrict pu8GmmPlusGray;
    u32 * restrict pu32MeanBg = (u32*)m_pfMeanBg;
    u32 u32MeanBg, u32Mask;
    CGMMPixel **pGmmTmp;
    l32 l32GmmPlusWidth, l32GmmPlusHeight;
    l32GmmPlusWidth  = m_l32Width / 2;
    l32GmmPlusHeight = m_l32Height / 2;
    pu8GmmPlusGray =  pu8YUV;
    
    //提取背景建模结果
    memset(pu8Mask, 0, m_l32Width * m_l32Height);
    pGmmTmp = m_ppcGmmpixel;
    u64Data = _itoll(0xFFFFFFFF, 0xFFFFFFFF);
    pu8OutBufferBase = pu8Mask;

    for(i = 0; i < m_l32NumBlocksH; i++)
    {
        for(j = 0; j < m_l32NumBlocksW; j++)
        {
            if(pGmmTmp[j]->fg == 255)           //foreground
            {
                pu8OutBuffer = &(pu8OutBufferBase[j * m_l32StepX]);
                
                //对输出8x8区域赋前景值
                for (m = 0; m < m_l32BlockHeight; m++)
                {
                    _mem8(pu8OutBuffer) = u64Data;
                    pu8OutBuffer += m_l32Width;
                }
            }
        }
        pu8OutBufferBase += m_l32StepY * m_l32Width; // pContext->iImgWidth;
        pGmmTmp += m_l32NumBlocksW;
    }
	//显示GMM结果
	ShowY8InOpenCV("GmmFg", pu8Mask, m_l32Width, m_l32Height, m_l32Width, 1);
    //---------------------GMM+------------------------
    
    //被定点化实现
    UpdateBgImg(pu8YUV, pu8MeanBg, pu8Mask);

     ShowY8InOpenCV("GmmBg: v4", pu8MeanBg, m_l32Width, m_l32Height, m_l32Width, 1);

    //求背景，前景的sobel边缘方向，然后降采样
    u8 * pu8VoteMapBase = (u8 *)m_ps16VData;
    int TgFix = m_f32Tg;
    int TsFix = m_f32Ts * 32;  
    int TdFix = 30;

    //高架设不易发生断裂，sobel边缘比较即可
    //if(TrackSettingInfo.m_HightMount) TdFix = 255;

	//原图和背景图的纹理和光照强度比较
    SobelCompare_n(pu8GmmPlusGray, pu8MeanBg, m_l32Width, m_l32Height, m_l32Width, 
        TgFix, TsFix, TdFix, pu8VoteMapBase, m_ps16HData, TRUE);

	//统计5x5块里面的投票结果
    Sum5x5(pu8VoteMapBase,  m_l32Width, m_l32Height, pu8VoteMapBase, (u8*)m_ps16HData);

	//如果大于12（5x5/2）则认为是前景。高架设取阈值为10
    ClearFalseFg(pu8Mask, pu8VoteMapBase, 12, m_l32Width, m_l32Height);

//static FILE * fp = NULL;
//if(fp==NULL) fp = fopen("J:\\test_video\\yuv\\lcgmmpp_5m.yuv","wb");
//if(fp) {fwrite(pu8Mask, 1, 704*396, fp);fflush(fp);}

    ShowY8InOpenCV("GmmPlusFg: v4", pu8Mask, m_l32Width, m_l32Height, m_l32Width, 2);

//     return TRUE;
}

void CLcgmm::process(u8 *pu8Img, u8 *pu8MeanBg, u8 *pu8Mask)
{
    //局部对比背景建模
    lcgmm(pu8Img); 

    //gmm+背景建模
    gmmplus(pu8Img, pu8MeanBg, pu8Mask);
}

//设置参数
void CLcgmm::setParameters(TLcgmmPlusParam *ptLcgmmPara)
{
    if(ptLcgmmPara)
    {
        m_tLcgmmPara.tStLcg = ptLcgmmPara->tStLcg;
        m_tLcgmmPara.tStGmplus = ptLcgmmPara->tStGmplus;

        m_f32Ts = ptLcgmmPara->tStGmplus.f32GmmPlusTs;      //0.4
        m_f32Tg = ptLcgmmPara->tStGmplus.f32GmmPlusTg * ptLcgmmPara->tStGmplus.f32GmmPlusTg;      //80x80
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

CLcgmm::CLcgmm(int l32Width, int l32Height):m_l32Width(l32Width),m_l32Height(l32Height)
{
    l32 l32BlocksH, l32BlocksW;
    l32 l32Index, l32Num;

    m_l32StepX = 4;
    m_l32StepY = 4;
    m_l32BlockWidth = 8;
    m_l32BlockHeight = 8;
    m_f32LearnRate = 0.005f; // Gmm+背景的学习率
    m_f32LcLearnRate = 0.01f; // 缩小图像的Gmm+背景的学习率

    m_l32NumBlocksW = int((l32Width - m_l32BlockWidth) / m_l32StepX);
    m_l32NumBlocksH = int((l32Height - m_l32BlockHeight) / m_l32StepY);

    l32BlocksH = m_l32NumBlocksH;
    l32BlocksW = m_l32NumBlocksW;

    m_ppcGmmpixel = new CGMMPixel*[l32BlocksW * l32BlocksH];
    m_ps8GmmBuf = (s8 *)AVMalloc(l32BlocksW * l32BlocksH * sizeof(CGMMPixel));

    if(NULL != m_ps8GmmBuf)
    {
        for(l32Index = 0; l32Index < l32BlocksW * l32BlocksH; l32Index++)
        {
            m_ppcGmmpixel[l32Index] = (CGMMPixel *)(m_ps8GmmBuf + sizeof(CGMMPixel) * l32Index);

            for(l32Num = 0; l32Num < K_MODELS; l32Num++)
            {
                m_ppcGmmpixel[l32Index]->model[l32Num].weight = 0.0f;
                m_ppcGmmpixel[l32Index]->model[l32Num].mean   = -999.0f;
                m_ppcGmmpixel[l32Index]->model[l32Num].var = InitialVariance;
                m_ppcGmmpixel[l32Index]->model[l32Num].dist2  = 0.0f;
                m_ppcGmmpixel[l32Index]->model[l32Num].mah2   = 0.0f; 
            }

            m_ppcGmmpixel[l32Index]->pmMatch = &m_ppcGmmpixel[l32Index]->model[0];
        }
    }

    m_ps16Iimg = (s16 *)AVMalloc((l32Width >> 2)*(l32Height >> 2)*sizeof(s16));
    m_pu32Iimg = (u32 *)AVMalloc((l32Width >> 2)*(l32Height >> 2)*sizeof(u32));

    // GMM+内存分配
    m_ps16VData = (s16 *)AVMalloc(l32Width * l32Height * sizeof(s16));
    m_ps16HData = (s16 *)AVMalloc(l32Width * l32Height * sizeof(s16));
    m_pfMeanBg = (float *)AVMalloc(l32Width * l32Height *  sizeof(float));

    if(NULL != m_pfMeanBg)
    {
        memset(m_pfMeanBg, 0, l32Width * l32Height * sizeof(float));
    }

    //背景建模初始化参数
    m_tLcgmmPara.tStGmplus.f32GmmPlusTs = 0.40000001;
    m_tLcgmmPara.tStGmplus.f32GmmPlusTg = 40;

    m_tLcgmmPara.tStLcg.f32AlphaPara = 0.002f;
    m_tLcgmmPara.tStLcg.f32RhoPara = 0.00375f;
    m_tLcgmmPara.tStLcg.f32SigmaPara = 2.5f;
    m_tLcgmmPara.tStLcg.f32ThreshPara = 0.12f;
    m_tLcgmmPara.tStLcg.f32PixelMaxVar = 0.01f;
    m_tLcgmmPara.tStLcg.f32PixelMinVar = 0.004f;	
	setParameters(&m_tLcgmmPara);
}

CLcgmm::~CLcgmm()
{
    if(NULL != m_ps16Iimg)
    {
        AVFree(m_ps16Iimg);
    }

    if(NULL != m_pu32Iimg)
    {
        AVFree(m_pu32Iimg);
    }

    if(NULL != m_ps8GmmBuf)
    {
        AVFree(m_ps8GmmBuf);
    }

    if(NULL != m_ppcGmmpixel)
    {
        delete [] m_ppcGmmpixel;
    }

    if(NULL != m_ps16VData)
    {
        AVFree(m_ps16VData);
    }

    if(NULL != m_ps16HData)
    {
        AVFree(m_ps16HData);
    }

    if(NULL != m_pfMeanBg)
    {
        AVFree(m_pfMeanBg);
    }
}

void CLcgmm::lcgmm(u8 *pu8YUV)
{
    int i,j,iStepY;
    float std;    
    float fLocalContrast;
    float myEx, myEx2;
    float tmp, tmp2;
    l32 l32IimgWidth;
    s16 *ps16IimgC, *ps16IimgP1;
    u32 *pu32ImgC, *pu32ImgP1;
    CGMMPixel **pGmmTmp;

    //-------------LCGMM start--------------------

    // 4cif灰度图像的4x4降采样局部求和
    iimg_4_n(pu8YUV, m_ps16Iimg, m_pu32Iimg, m_l32Height, m_l32Width);

    //循环中的局部变量和局部指针初始化
    tmp = _rcpsp(64.0);
    pGmmTmp = m_ppcGmmpixel;
    l32IimgWidth = m_l32Width >> 2;
    pu32ImgC = m_pu32Iimg;
    pu32ImgP1 = m_pu32Iimg + l32IimgWidth;
    ps16IimgC = m_ps16Iimg;
    ps16IimgP1 = m_ps16Iimg + l32IimgWidth;

    for(i = 0; i < m_l32NumBlocksH; i++)
    {
        iStepY = i*m_l32StepY;
        for(j = 0; j < m_l32NumBlocksW; j++)
        {
            //compute local contrast
            myEx = (ps16IimgC[j] + ps16IimgC[j + 1] + ps16IimgP1[j] + ps16IimgP1[j + 1]) * tmp;
            myEx2 = (pu32ImgC[j] + pu32ImgC[j + 1] + pu32ImgP1[j] + pu32ImgP1[j + 1]) * tmp;

            std = _rsqrsp(myEx2  - myEx * myEx);
            std = _rcpsp(std);

            if (std<2.0)
            {
                std=2.0;
            }
            tmp2 = _rcpsp(myEx+5.0);
            fLocalContrast=(std * tmp2);

            //GMM
            pGmmTmp[j]->Process2(fLocalContrast, m_tLcgmmPara.tStLcg);

            //人为去除对比度过低处前景检测的结果
            if(!(fLocalContrast>m_tLcgmmPara.tStLcg.f32ThreshPara))
            {
                pGmmTmp[j]->fg = 0;
            }
        }

        ps16IimgC += l32IimgWidth;
        ps16IimgP1 += l32IimgWidth;
        pu32ImgC += l32IimgWidth;
        pu32ImgP1 += l32IimgWidth;
        pGmmTmp += m_l32NumBlocksW;
    }
}

void CLcgmm::UpdateBgImg(u8 *pu8YUV, u8 *pu8MeanBg, u8 *pu8Mask) // VPM_IMGPROCESS_CONTEXT_STRUCT *pContext)
{
    update_bg_img(pu8YUV, pu8MeanBg, pu8Mask, (u16 *)m_pfMeanBg, m_l32Width, m_l32Height, m_f32LcLearnRate);
}



///////////////////////////////////////////////////////////////////////////////////////
CGMMPixel::CGMMPixel()
{
    for(int i=0;i<K_MODELS;i++)
    {
        model[i].weight = 0.0f;
        model[i].mean   = -999.0f;
        model[i].var = InitialVariance;
        model[i].dist2  = 0.0f;
        model[i].mah2   = 0.0f; 
    }

    pmMatch = &model[0];     
}

CGMMPixel::~CGMMPixel()
{
}

const int CompareTable[8][3]=
{
    {0,1,2},//
    {0,2,1},
    {1,0,2},
    {0,1,2},
    {2,1,0},
    {2,0,1},
    {1,2,0},
    {0,1,2}//
};

void CGMMPixel::Process2(float data, TLcgmmParam &Param)
{
    int i, j, compid, match_i;
    PixelModel pm;
    float tmp, sumw;
    //---------------------sort-------------------------
    //---------------------start-------------------------
    float fsigma2 = Param.f32SigmaPara * Param.f32SigmaPara;
#if 1
    float v0, v1, v2;
    tmp = _rcpsp(model[0].var);
    v0 = model[0].weight*model[0].weight * tmp;

    tmp = _rcpsp(model[1].var);
    v1 = model[1].weight*model[1].weight * tmp;

    tmp = _rcpsp(model[2].var);
    v2 = model[2].weight*model[2].weight * tmp ;
    compid = (v0 > v1) + (v1 > v2)*2 + (v2 > v0)*4;

    // see if the current pixel matches any of the models
    sumw = 0;
    for(i=0; i<K_MODELS; i++)
    {
        match_i = CompareTable[compid][i];

        model[match_i].dist2 = (data - model[match_i].mean)*(data - model[match_i].mean);
        tmp = _rcpsp(model[match_i].var);
        model[match_i].mah2 = model[match_i].dist2 * tmp; 
        if (model[match_i].mah2 < fsigma2) 
            break;

        sumw += model[match_i].weight;
    }
#else
    float v[K_MODELS];
    for(i=0; i<K_MODELS; i++)
    {
        v[i] = model[i].weight*model[i].weight / model[i].var;    
    }
    for(i=0; i<K_MODELS-1; i++)
    {
        for(j=i+1; j<K_MODELS; j++)
        {
            if (v[i] < v[j])
            {
                tmp = v[i]; v[i] = v[j]; v[j] = tmp;
                pm = model[i]; 
                model[i] = model[j]; 
                model[j] = pm;
            }
        }
    }
    //------------------------end------------------------
    // see if the current pixel matches any of the models
    sumw = 0;
    for(i=0; i<K_MODELS; i++)
    {
        model[i].dist2 = (data - model[i].mean)*(data - model[i].mean);
        model[i].mah2 = model[i].dist2 / model[i].var; 
        if (model[i].mah2 < sigma*sigma) 
            break;

        sumw += model[i].weight;
    }
#endif
    if (i<K_MODELS)
    {
        // we matched an existing model
        pmMatch = &model[match_i];
        //---------------------end-------------------------

        // find which models account for the bg by summing the model's weights until we exceed <bgT>

#if 1
        //判断命中背景部分还是前景部分
        if(sumw > bgT)  fg=255;
        else fg=0;

        // now that we've classified, we need to update the models
        // zot: do we also update the weights if nothing matched?
        // update the weights (we normalize later)
        tmp = 1.0f - Param.f32AlphaPara;

        for(i=0; i<K_MODELS; i++)
        {
            model[i].weight *= tmp;
        }
        pmMatch->weight += Param.f32AlphaPara;
        for(i=0; i<K_MODELS; i++)
        {
            //解决贞子效应
            if (model[i].weight < 0.00001f)
                model[i].weight=0.0;
        }
#else

        tmp = 0.0;    //store sum of weights 
        for(j=0; j<K_MODELS; j++)
        {
            tmp += model[j].weight;
            if (tmp > bgT) break;
        }
        // determine if the matched model is part of the bg or fg
        fg=0;
        for(i=j+1; i<K_MODELS; i++)
        {
            if (pmMatch == &model[i])
            {
                fg=255;
                break;
            }
        }
        // now that we've classified, we need to update the models
        // zot: do we also update the weights if nothing matched?
        // update the weights (we normalize later)
        tmp = 1.0f - alpha;
        for(i=0; i<K_MODELS; i++)
        {
            if (&model[i] == pmMatch)
            {
                // this is the matched model, so M = 1
                model[i].weight = tmp * model[i].weight + alpha;
            }
            else{
                // this is not the matched model, so M = 0
                model[i].weight = tmp * model[i].weight;
            }
            //解决贞子效应
            if (model[i].weight < 0.00001f)
                model[i].weight=0.0;

        }

#endif

        // update the parameters for the matched model
        //rho = alpha * Gauss(pmMatch);
        tmp = 1.0f - Param.f32RhoPara;

        // first we update the mean
        pmMatch->mean = (tmp * pmMatch->mean) + (Param.f32RhoPara * data);
        // and then we update the variance
        pmMatch->var  = (tmp * pmMatch->var) + (Param.f32RhoPara * pmMatch->dist2);

        // not to make var too small or too large
        if(pmMatch->var < Param.f32PixelMinVar)
            pmMatch->var = Param.f32PixelMinVar; 
        if(pmMatch->var > Param.f32PixelMaxVar)
            pmMatch->var = Param.f32PixelMaxVar; 
    }
    else{
        // we didn't match anything 
        //--------------------------start----------------------
        tmp = model[0].weight;
        j=0;
        for(i=1; i<K_MODELS; i++)
        {
            if (model[i].weight <= tmp)
            {
                tmp = model[i].weight;
                j = i;
            }
        }
        pmMatch = &model[j];
        //-------------------------end--------------------------

        // replace least probable model with the current data
        pmMatch->weight = InitialWeight;
        pmMatch->mean   = data;
        pmMatch->var    = InitialVariance;
        // we assume the 'new' model is foreground
        fg=255;
    }   
    // we adjusted the weights during the update, so we need to normalize
    //--------------------------start-----------------------------
    tmp=0.0;
    tmp = model[0].weight + model[1].weight + model[2].weight;

    tmp = _rcpsp(tmp);

    model[0].weight *= tmp;
    model[1].weight *= tmp;
    model[2].weight *= tmp;

    //--------------------------end---------------------------

}
//////////////////////////////////////////////////////////////////////
void CGMMPixel::Process(float data, TLcgmmParam &Param)
{
    int i, j;
    float v[3];
    PixelModel pm;
    float tmp;

    //---------------------sort-------------------------
    //---------------------start-------------------------
    for(i=0; i<K_MODELS; i++)
    {
        v[i] = model[i].weight*model[i].weight / model[i].var;    
    }

    for(i=0; i<K_MODELS-1; i++)
    {
        for(j=i+1; j<K_MODELS; j++)
        {
            if (v[i] < v[j])
            {
                tmp = v[i]; v[i] = v[j]; v[j] = tmp;
                pm = model[i]; 
                model[i] = model[j]; 
                model[j] = pm;
            }
        }
    }
    //------------------------end------------------------
    // see if the current pixel matches any of the models

    for(i=0; i<K_MODELS; i++)
    {
        model[i].dist2 = (data - model[i].mean)*(data - model[i].mean);
        model[i].mah2 = model[i].dist2 / model[i].var; 
        if (model[i].mah2 < Param.f32SigmaPara * Param.f32SigmaPara) 
        {
            break;
        }
    }

    if (i<K_MODELS)
    {
        // we matched an existing model
        pmMatch = &model[i];
        //---------------------end-------------------------
        // find which models account for the bg by summing the model's weights until we exceed <bgT>
        tmp = 0.0;    //store sum of weights 
        for(j=0; j<K_MODELS; j++)
        {
            tmp += model[j].weight;
            if (tmp > bgT) break;
        }

        // determine if the matched model is part of the bg or fg
        fg=0;
        for(i=j+1; i<K_MODELS; i++)
        {
            if (pmMatch == &model[i])
            {
                fg=255;
                break;
            }
        }

        // now that we've classified, we need to update the models
        // zot: do we also update the weights if nothing matched?
        // update the weights (we normalize later)
        tmp = 1.0f - Param.f32AlphaPara;
        for(i=0; i<K_MODELS; i++)
        {
            if (&model[i] == pmMatch)
            {
                // this is the matched model, so M = 1
                model[i].weight = tmp * model[i].weight + Param.f32AlphaPara;
            }
            else{
                // this is not the matched model, so M = 0
                model[i].weight = tmp * model[i].weight;
            }
            //解决贞子效应
            if (model[i].weight < 0.00001f)
            {
                model[i].weight=0.0;
            }
        }

        // update the parameters for the matched model
        //rho = alpha * Gauss(pmMatch);
        tmp = 1.0f - Param.f32RhoPara;

        // first we update the mean
        pmMatch->mean = (tmp * pmMatch->mean) + (Param.f32RhoPara * data);
        // and then we update the variance
        pmMatch->var  = (tmp * pmMatch->var) + (Param.f32RhoPara * pmMatch->dist2);

        // not to make var too small or too large
        if(pmMatch->var < Param.f32PixelMinVar)
            pmMatch->var = Param.f32PixelMinVar; 
        if(pmMatch->var > Param.f32PixelMaxVar)
            pmMatch->var = Param.f32PixelMaxVar; 
    }
    else
    {
        // we didn't match anything 
        //--------------------------start----------------------
        tmp = model[0].weight;
        j=0;
        for(i=1; i<K_MODELS; i++)
        {
            if (model[i].weight <= tmp)
            {
                tmp = model[i].weight;
                j = i;
            }
        }
        pmMatch = &model[j];
        //-------------------------end--------------------------

        // replace least probable model with the current data
        pmMatch->weight = InitialWeight;
        pmMatch->mean   = data;
        pmMatch->var    = InitialVariance;
        // we assume the 'new' model is foreground
        fg=255;
    }   
    // we adjusted the weights during the update, so we need to normalize
    //--------------------------start-----------------------------
    tmp=0.0;
    for(i=0; i<K_MODELS; i++) 
    {
        tmp += model[i].weight;
    }

    for(i=0; i<K_MODELS; i++) 
    {
        model[i].weight /= tmp;
    }
    //--------------------------end---------------------------
}
