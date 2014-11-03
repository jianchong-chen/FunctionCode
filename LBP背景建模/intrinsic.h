/*****************************************************************************
  模块名      : intrinsic.h
  文件名      : intrinsic.h
  相关文件    : 
  文件实现功能: VC下模拟TI intrinsic的内联函数实现
  作者        : 
  版本        : V1.0
******************************************************************************/

#ifndef _INTRINSIC_H_
#define _INTRINSIC_H_

#ifdef ARCH_DM64XPLUS

#define BRKAT(a)
#define PROFILER_DEF(a)
#define PROFILER_CALL(a, c)

static INLINE u32 Bswap32(u32 u32X)
{
    u32X = _swap4(u32X);
    u32X = _packlh2(u32X, u32X);
    return u32X;
}

//CCS下右移由于使用的是40bit的移位,因此右移32位也能得到合法的
//结果0
#define _SHRU(a,b)  ((a) >> (b))

typedef unsigned long u40;//40bit integer

#else

#pragma warning( disable : 4142 )
#pragma warning( disable : 4028 )

//////////////////////////////////////////////////////////////////////////
//下面的函数仅供VC中模拟统计TI DSP代码调用次数,开销使用.
typedef struct _CYCLE_PROFILER_
{
    int total_cycle;
    int call_times;
}CYCLE_PROFILER;

#define PROFILER_DEF(a) CYCLE_PROFILER a={0,0};
#define PROFILER_CALL(a, c) do{a.call_times++; a.total_cycle+=c;}while(0);
//////////////////////////////////////////////////////////////////////////


#define BYTE0(u32A) ((u32)(u32A) & 0xFF)
#define BYTE1(u32A) ((u32)(u32A >> 8) & 0xFF)
#define BYTE2(u32A) ((u32)(u32A >> 16) & 0xFF)
#define BYTE3(u32A) ((u32)(u32A) >> 24)

#define SBYTE0(u32A) (((l32)(u32A<<24)) >> 24)
#define SBYTE1(u32A) (((l32)(u32A<<16)) >> 24)
#define SBYTE2(u32A) (((l32)(u32A<<8)) >> 24)
#define SBYTE3(u32A) (((l32)(u32A)) >> 24)

#define WORD0(u32A) ((u32)(u32A) & 0xFFFF)
#define WORD1(u32A) ((u32)(u32A) >> 16)

#define BRKAT(a) do{if(a) __asm{int 3}; }while(0);

#define restrict

#define T_ASSERT(a)  while(!(a)){}
#include <math.h>
static __inline u32 Bswap32(u32 u32X)
{
    u32X = ((u32X << 8) & 0xFF00FF00) | ((u32X >> 8) & 0x00FF00FF);
    return (u32X >> 16) | (u32X << 16);
}

//这个函数模仿TI DSP的LMBD 1 的情况，返回从左向右找到的第一个比特1
//的位置，
#define _my_lmbd1 _lmbd1

static __inline s64 _itoll(u32 u32Src1, u32 u32Src2)
{
    u64 u64Res;

    u64Res = u32Src1;
    u64Res = (u64Res << 32) | u32Src2;

    return u64Res;
}

static __inline l32 _lmbd1(u32 u32Reg)
{
	l32 l32Index;
	for(l32Index = 31; l32Index >= 0; l32Index--)
	{
		if(u32Reg & (1 << l32Index))
		{
			return (31 - l32Index);
		}
	}
	//未找到，返回32
	return 32;
}

static __inline l32 _lmbd0(u32 u32Reg)
{
	l32 l32Index;
	for(l32Index = 31; l32Index >= 0; l32Index--)
	{
		if((u32Reg & (1 << l32Index)) == 0)
		{
			return (31 - l32Index);
		}
	}
	//未找到，返回32
	return 32;
}
#if ((defined(ARCH_X86_WIN32)) ||  (defined(ARCH_X86_LINUX)) || (defined(ARCH_ATOM_WIN32)) || (defined(ARCH_ATOM_LINUX)))
static __inline int _lmbd(l32 l32BitDetect,u32 u32Reg)
{
    l32 l32Index;
	if(l32BitDetect == 1) 
    {
        //使用x86对应指令会更有效一点
        l32Index = -1;
        __asm
        {
            bsr eax, u32Reg 
            jz  _u32RegIsZero
            mov l32Index, eax
_u32RegIsZero:
        };
        
        return 31 - l32Index;

        //return _lmbd1(u32Reg);
    }
	else 
    {
        return _lmbd0(u32Reg);
    }
}
#else
static __inline int _lmbd(l32 l32BitDetect,u32 u32Reg)
{
	if(l32BitDetect == 1) 
	{
		return _lmbd1(u32Reg);
	}
	else 
	{
		return _lmbd0(u32Reg);
	}
}
#endif

static __inline  u32 _extu(u32 u32Src, u32 u32Csta, u32 u32Cstb)
{
	T_ASSERT(u32Cstb<32);
	return (u32Src << u32Csta >> u32Cstb);
}

static __inline  l32 _ext(l32 l32Src, u32 u32Csta, u32 u32Cstb)
{
	T_ASSERT(u32Cstb<32);
	return (l32Src << u32Csta >> u32Cstb);
}
static __inline  u32 _max2(l32 l32a, l32 l32b)
{
	s16 s16Low;
	u32 u32Max;

	if((s16)(l32a & 0xFFFF) > (s16)(l32b & 0xFFFF))
	{
		s16Low = (s16)(l32a & 0xFFFF);
	}
	else
	{
		s16Low = (s16)(l32b & 0xFFFF);
	}

	if((l32a >> 16) > (l32b >> 16) )
	{
		u32Max = (l32a >> 16);
	}
	else
	{
		u32Max = (l32b >> 16);
	}
	return (u32Max << 16) | (u32)(s16Low & 0xFFFF);
}
static __inline  u32 _min2(l32 l32a, l32 l32b)
{
	s16 s16Low;
	u32 u32Max;

	if((s16)(l32a & 0xFFFF) < (s16)(l32b & 0xFFFF))
	{
		s16Low = (s16)(l32a & 0xFFFF);
	}
	else
	{
		s16Low = (s16)(l32b & 0xFFFF);
	}

	if((l32a >> 16) < (l32b >> 16) )
	{
		u32Max = (l32a >> 16);
	}
	else
	{
		u32Max = (l32b >> 16);
	}
	return (u32Max << 16) | (u32)(s16Low & 0xFFFF);
}

static __inline  u32 _pack2(u32 u32Src1, u32 u32Src2)
{
    return (u32Src1 << 16) | (u32Src2 & 0x0000FFFF);
}
static __inline  u32 _packh2(u32 u32Src1, u32 u32Src2)
{
    return (u32Src1 & 0xFFFF0000) | (u32Src2 >> 16);
}
static __inline  u32 _packl4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret;
    
    u32Ret = (u32Src1 >> 16) & 0xFF;
    u32Ret = (u32Ret << 8) | (u32Src1 & 0xFF);
    u32Ret = (u32Ret << 8) | ((u32Src2 >> 16) & 0xFF);
    u32Ret = (u32Ret << 8) | (u32Src2 & 0xFF);
    
    return u32Ret;
}
static __inline  u32 _packh4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret;
    
    u32Ret = (u32Src1 >> 24) & 0xFF;
    u32Ret = (u32Ret << 8) | ((u32Src1 >> 8) & 0xFF);
    u32Ret = (u32Ret << 8) | ((u32Src2 >> 24) & 0xFF);
    u32Ret = (u32Ret << 8) | ((u32Src2 >> 8) & 0xFF);
    
    return u32Ret;
}

static __inline  u32 _clr(u32 u32Src, u32 u32Csta, u32 u32Cstb)
{
    u32 u32MaskSize;
    u32 u32Mask;
    u32MaskSize = (u32Cstb - u32Csta + 1);
    if(u32MaskSize == 32)
    {
        u32Mask = 0xFFFFFFFF;
    }
    else
    {
        u32Mask = ((u32)1 << u32MaskSize) - 1;
        u32Mask <<= u32Csta;
    }
    return (u32Src & (~u32Mask));
}

static __inline  u32 _shru_ti(u32 u32Src, u32 u32ShiftAmount)
{
	if(u32ShiftAmount > 31)
	{
		return 0;
	}
	else
	{
		return (u32Src >> u32ShiftAmount);
	}
}

static __inline  x64 _itod(u32 u32hi, u32 u32lo)
{
	u64 x64Ret;
	x64Ret = u32hi;
	x64Ret = (x64Ret << 32) | u32lo;
	return x64Ret;
}

static __inline  u32 _bitr(u32 u32Src)
{
    u32 u32Dst = 0;
    l32 l32Index;
    for(l32Index = 0 ; l32Index < 32 ; l32Index ++)
    {
        if(u32Src & (1 << l32Index))
        {
            u32Dst |= ((u32)(0x80000000) >> l32Index);
        }
    }
    return u32Dst;
}
static __inline  l32 _dotp2(l32 l32A, l32 l32B)
{
    return _ext(l32A, 0, 16) * _ext(l32B, 0, 16) +
           _ext(l32A, 16, 16) * _ext(l32B, 16, 16); 
}
static __inline  l32 _dotpn2(l32 l32A, l32 l32B)
{
    return _ext(l32A, 0, 16) * _ext(l32B, 0, 16) - 
           _ext(l32A, 16, 16) * _ext(l32B, 16, 16); 
}
static __inline  l32 _dotprsu2(l32 l32A, u32 u32B)
{
    return (_ext(l32A, 0, 16) * _extu(u32B, 0, 16) +
            _ext(l32A, 16, 16) * _extu(u32B, 16, 16) + 0x8000) >> 16; 
}

static __inline  x64 _addsub2(l32 l32Src32, l32 l32Src10)
{
    u32 u32AddDst32,u32AddDst10;
    u32 u32SubDst32,u32SubDst10;

    u32AddDst10 = _ext(l32Src32, 16, 16) + _ext(l32Src10, 16, 16);
    u32AddDst32 = _ext(l32Src32, 0, 16) + _ext(l32Src10, 0, 16);
    
    u32SubDst10 = _ext(l32Src32, 16, 16) - _ext(l32Src10, 16, 16);
    u32SubDst32 = _ext(l32Src32, 0, 16) - _ext(l32Src10, 0, 16);

    return _itod( _pack2(u32AddDst32, u32AddDst10),_pack2(u32SubDst32, u32SubDst10) );
}

static __inline  u32 _cmpeq4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret = 0;

    if((u32Src1 & 0xFF)==(u32Src2 & 0xFF))
    {
        u32Ret |= 1;
    }
    if((u32Src1 & 0xFF00)==(u32Src2 & 0xFF00))
    {
        u32Ret |= 2;
    }
    if((u32Src1 & 0xFF0000)==(u32Src2 & 0xFF0000))
    {
        u32Ret |= 4;
    }
    if((u32Src1 & 0xFF000000)==(u32Src2 & 0xFF000000))
    {
        u32Ret |= 8;
    }
    return u32Ret;
}

static __inline  u32 _dcmpeq4(u64 src1, u64 src2)
{
    u32 u32Ret1 = _cmpeq4((u32)src1, (u32)src2);
    u32 u32Ret2 = _cmpeq4((u32)(src1>>32), (u32)(src2>>32));
    return u32Ret1 | u32Ret2 << 4;
}

static __inline  u32 _cmpgtu4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret = 0;
    
    if((u32Src1 & 0xFF) > (u32Src2 & 0xFF))
    {
        u32Ret |= 1;
    }
    if((u32Src1 & 0xFF00) > (u32Src2 & 0xFF00))
    {
        u32Ret |= 2;
    }
    if((u32Src1 & 0xFF0000) > (u32Src2 & 0xFF0000))
    {
        u32Ret |= 4;
    }
    if((u32Src1 & 0xFF000000) > (u32Src2 & 0xFF000000))
    {
        u32Ret |= 8;
    }
    return u32Ret;
}
static __inline  u32 _satu8(l32 l32Temp)
{
    if(l32Temp < 0) 
    {
        l32Temp = 0;
    }
    if(l32Temp > 255) 
    {
        l32Temp = 255;
    }
    return l32Temp;
}
static __inline  l32 _sats16(l32 l32Temp)
{
    if(l32Temp < ((-1)<<15)) 
    {
        l32Temp = ((-1)<<15);
    }
    if(l32Temp > (0x7FFF)) 
    {
        l32Temp = 0x7FFF;
    }
    return l32Temp;
}
static __inline  u32 _spack2(l32 l32Src1, l32 l32Src2)
{
    u32 u32Ret = 0;
    u32Ret = _sats16(l32Src1) << 16;
    u32Ret |= _sats16(l32Src2);
    return u32Ret;
}
static __inline  u32 _spacku4(l32 l32Src1, l32 l32Src2)
{
    u32 u32Ret = 0;
    l32 l32Temp;

    l32Temp = _ext(l32Src1, 0, 16);
    u32Ret |= _satu8(l32Temp) << 24;

    l32Temp = _ext(l32Src1, 16, 16);
    u32Ret |= _satu8(l32Temp) << 16;

    l32Temp = _ext(l32Src2, 0, 16);
    u32Ret |= _satu8(l32Temp) << 8;
    
    l32Temp = _ext(l32Src2, 16, 16);
    u32Ret |= _satu8(l32Temp);
    
    return u32Ret;
}
static __inline  u32 _shfl(u32 u32Src1)
{
    l32 l32Index;
    u32 u32Res = 0;
    u32 u32D1 = u32Src1 >> 16;
    u32 u32D0 = u32Src1 & 0xFFFF;
    for (l32Index = 0; l32Index < 16; l32Index++)
    {
        u32Res |= ((u32D0 >> l32Index) & 1) << (l32Index*2 + 0);
        u32Res |= ((u32D1 >> l32Index) & 1) << (l32Index*2 + 1);
    }
    return u32Res;
}
static __inline  x64 _shfl3(l32 l32Src1, l32 l32Src2)
{
    l32 l32Index;
    l32 l32Inp0, l32Inp1, l32Inp2;
    u64 x64Result;
    l32Inp0 = l32Src2 & 0xFFFF;
    l32Inp1 = l32Src1 & 0xFFFF;
    l32Inp2 = (l32Src1 >> 16) & 0xFFFF;
    x64Result = 0;
    for (l32Index = 0; l32Index < 16; l32Index++)
    {
        x64Result |=  (u64)(((l32Inp0 >> l32Index) & 1)) << ((l32Index * 3) + 0);
        x64Result |=  (u64)(((l32Inp1 >> l32Index) & 1)) << ((l32Index * 3) + 1);
        x64Result |=  (u64)(((l32Inp2 >> l32Index) & 1)) << ((l32Index * 3) + 2);
    }
    return x64Result;
}
static __inline u32 _packlh2(l32 l32Src1, l32 l32Src2)
{
    u32 u32Ret;
    u32Ret =  (l32Src1 << 16);
    u32Ret |= (_extu(l32Src2, 0, 16));
    return u32Ret;
}
static __inline u32 _packhl2(l32 l32Src1, l32 l32Src2)
{
    u32 u32Ret;
    u32Ret =  (l32Src1 & 0xFFFF0000);
    u32Ret |= (l32Src2 & 0x0000FFFF);
    return u32Ret;
}
static __inline  u32 _add2(l32 l32Src1, l32 l32Src2 )
{
    u32 u32Ret;
    u32Ret =  ((_ext(l32Src1, 0, 16) + _ext(l32Src2, 0, 16)) & 0xFFFF) << 16;
    u32Ret |= (_ext(l32Src1, 16, 16) + _ext(l32Src2, 16, 16)) & 0xFFFF;
    return u32Ret;    
}
static __inline  u32 _sadd2(l32 l32Src1, l32 l32Src2 )
{
    u32 u32Ret;
    l32 l32Tmp;

    l32Tmp = _ext(l32Src1, 0, 16) + _ext(l32Src2, 0, 16);
    if(l32Tmp < -32768) l32Tmp = -32768;
    if(l32Tmp > 32767) l32Tmp = 32767;
    u32Ret  =  (l32Tmp) << 16;

    l32Tmp = _ext(l32Src1, 16, 16) + _ext(l32Src2, 16, 16);
    if(l32Tmp < -32768) l32Tmp = -32768;
    if(l32Tmp > 32767) l32Tmp = 32767;
    u32Ret |= (l32Tmp) & 0xFFFF;
    return u32Ret;      
}

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : (-(x)))
#endif

static __inline  u32 _abs2(l32 l32Src1)
{
    u32 u32Ret;
    u32Ret =  ABS(_ext(l32Src1, 0, 16)) << 16;
    u32Ret |= ABS(_ext(l32Src1, 16, 16));
    return u32Ret;
}

static __inline l32 _subc(l32 l32Src1, l32 l32Src2)
{
    l32 l32Ret;
    if(l32Src1 >= l32Src2)
    {
        l32Ret = ((l32Src1-l32Src2)<<1) + 1;
    }
    else
    {
        l32Ret = l32Src1 << 1;
    }
    return l32Ret;
}

static __inline  u32 _xpnd2(l32 l32Src1)
{
    u32 u32Ret = 0;
    if(l32Src1 & 1)
    {
        u32Ret |= 0xFFFF;
    }
    if(l32Src1 & 2)
    {
        u32Ret |= 0xFFFF0000;
    }
    return u32Ret;
}

static __inline  u32 _xpnd4(l32 l32Src1)
{
    u32 u32Ret = 0;
    if(l32Src1 & 1)
    {
        u32Ret |= 0xFF;
    }
    if(l32Src1 & 2)
    {
        u32Ret |= 0xFF00;
    }
    if(l32Src1 & 4)
    {
        u32Ret |= 0xFF0000;
    }
    if(l32Src1 & 8)
    {
        u32Ret |= 0xFF000000;
    }
    return u32Ret;
}

static __inline x64 _mpysu4(u32 u32Src1, u32 u32Src2)
{
    u64 x64ret=0;
    
    x64ret |= ((SBYTE3(u32Src1) * BYTE3(u32Src2)) & 0xFFFF);
	x64ret <<= 16;
    x64ret |= ((SBYTE2(u32Src1) * BYTE2(u32Src2)) & 0xFFFF);
	x64ret <<= 16;
    x64ret |= ((SBYTE1(u32Src1) * BYTE1(u32Src2)) & 0xFFFF);
	x64ret <<= 16;
    x64ret |= ((SBYTE0(u32Src1) * BYTE0(u32Src2)) & 0xFFFF);
    
    return x64ret;
}

static __inline x64 _mpyu4(u32 u32Src1, u32 u32Src2)
{
    u64 x64ret=0;
    
    x64ret |= ((BYTE3(u32Src1) * BYTE3(u32Src2)) & 0xFFFF);
	x64ret <<= 16;
    x64ret |= ((BYTE2(u32Src1) * BYTE2(u32Src2)) & 0xFFFF);
	x64ret <<= 16;
    x64ret |= ((BYTE1(u32Src1) * BYTE1(u32Src2)) & 0xFFFF);
	x64ret <<= 16;
    x64ret |= ((BYTE0(u32Src1) * BYTE0(u32Src2)) & 0xFFFF);
    
    return x64ret;
}

static __inline  u32 _cmpgt2(l32 l32Src1, l32 l32Src2)
{
    u32 u32Ret = 0;
    
    if(_ext(l32Src1, 16, 16)>_ext(l32Src2, 16, 16))
    {
        u32Ret |= 1;
    }
    if(_ext(l32Src1, 0, 16)>_ext(l32Src2, 0, 16))
    {
        u32Ret |= 2;
    }
    return u32Ret;
}
static __inline  u32 _cmpeq2(l32 l32Src1, l32 l32Src2)
{
    u32 u32ret = 0;
    
    if(_ext(l32Src1, 16, 16) == _ext(l32Src2, 16, 16))
    {
        u32ret |= 1;
    }
    if(_ext(l32Src1, 0, 16) == _ext(l32Src2, 0, 16))
    {
        u32ret |= 2;
    }
    return u32ret;
}

static __inline  u32 _sub2(u32 u32Src1, u32 u32Src2 )
{
    u32 u32Ret;
    u32Ret  =  ((_ext(u32Src1, 0, 16) - _ext(u32Src2, 0, 16)) & 0xFFFF) << 16;
    u32Ret |= (_ext(u32Src1, 16, 16) - _ext(u32Src2, 16, 16)) & 0xFFFF;
    return u32Ret;    
}
static __inline  u32 _ssub2(u32 u32Src1, u32 u32Src2 )
{
    u32 u32Ret;
    l32 l32Tmp;

    l32Tmp = _ext(u32Src1, 0, 16) - _ext(u32Src2, 0, 16);
    if(l32Tmp < -32768) l32Tmp = -32768;
    if(l32Tmp > 32767) l32Tmp = 32767;
    u32Ret  =  (l32Tmp) << 16;

    l32Tmp = _ext(u32Src1, 16, 16) - _ext(u32Src2, 16, 16);
    if(l32Tmp < -32768) l32Tmp = -32768;
    if(l32Tmp > 32767) l32Tmp = 32767;
    u32Ret |= (l32Tmp) & 0xFFFF;
    return u32Ret;    
}
static __inline  u32 _shr2(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret;
    u32Ret =  ((_ext(u32Src1, 0, 16) >> u32Src2) & 0xFFFF) << 16;
    u32Ret |= (_ext(u32Src1, 16, 16) >> u32Src2) & 0xFFFF;
    return u32Ret;
}
static __inline  u32 _shru2(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret;
    u32Ret =  ((_extu(u32Src1, 0, 16) >> u32Src2) & 0xFFFF) << 16;
    u32Ret |= (_extu(u32Src1, 16, 16) >> u32Src2) & 0xFFFF;
    return u32Ret;
}

static __inline  u32 _unpklu4(u32 u32Src1)
{
    u32 u32Ret = 0;
    u32Ret |= _extu(u32Src1, 16, 24) << 16;
    u32Ret |= _extu(u32Src1, 24, 24);
    return u32Ret;
}
static __inline  u32 _unpkhu4(u32 u32Src1)
{
    u32 u32Ret = 0;
    u32Ret |= _extu(u32Src1, 0, 24) << 16;
    u32Ret |= _extu(u32Src1, 8, 24);
    return u32Ret;
}

static __inline  u32 _dotpu4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret = 0;
    
    u32Ret += BYTE0(u32Src1) * BYTE0(u32Src2);
    u32Ret += BYTE1(u32Src1) * BYTE1(u32Src2);
    u32Ret += BYTE2(u32Src1) * BYTE2(u32Src2);
    u32Ret += BYTE3(u32Src1) * BYTE3(u32Src2);
    
    return u32Ret;
}
static __inline  l32 _dotpsu4(u32 u32Src1, u32 u32Src2)
{
    l32 l32Ret = 0;
    
    l32Ret += SBYTE0(u32Src1) * BYTE0(u32Src2);
    l32Ret += SBYTE1(u32Src1) * BYTE1(u32Src2);
    l32Ret += SBYTE2(u32Src1) * BYTE2(u32Src2);
    l32Ret += SBYTE3(u32Src1) * BYTE3(u32Src2);
    
    return l32Ret;
}

static __inline  u32 _subabs4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Ret = 0;
    l32 l32Diff0;
    l32 l32Diff1;
    l32 l32Diff2;
    l32 l32Diff3;
    
    l32Diff0 = (l32)BYTE0(u32Src1) - (l32)BYTE0(u32Src2);
    l32Diff1 = (l32)BYTE1(u32Src1) - (l32)BYTE1(u32Src2);
    l32Diff2 = (l32)BYTE2(u32Src1) - (l32)BYTE2(u32Src2);
    l32Diff3 = (l32)BYTE3(u32Src1) - (l32)BYTE3(u32Src2);

    u32Ret |= (ABS(l32Diff0) & 0xFF);
    u32Ret |= (ABS(l32Diff1) & 0xFF) << 8;
    u32Ret |= (ABS(l32Diff2) & 0xFF) << 16;
    u32Ret |= (ABS(l32Diff3) & 0xFF) << 24;
    return u32Ret;
}
static __inline  u32 _maxu4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Max0;
    u32 u32Max1;
    u32 u32Max2;
    u32 u32Max3;

    u32Max0 = (u32Src1 & 0xFF) > (u32Src2 & 0xFF)? (u32Src1 & 0xFF):(u32Src2 & 0xFF);
    u32Max1 = (u32Src1 & 0xFF00) > (u32Src2 & 0xFF00)? (u32Src1 & 0xFF00):(u32Src2 & 0xFF00);
    u32Max2 = (u32Src1 & 0xFF0000) > (u32Src2 & 0xFF0000)? (u32Src1 & 0xFF0000):(u32Src2 & 0xFF0000);
    u32Max3 = (u32Src1 & 0xFF000000) > (u32Src2 & 0xFF000000)? (u32Src1 & 0xFF000000):(u32Src2 & 0xFF000000);

    return u32Max0|u32Max1|u32Max2|u32Max3;
}
static __inline  u32 _minu4(u32 u32Src1, u32 u32Src2)
{
    u32 u32Min0;
    u32 u32Min1;
    u32 u32Min2;
    u32 u32Min3;

    u32Min0 = (u32Src1 & 0xFF) < (u32Src2 & 0xFF)? (u32Src1 & 0xFF):(u32Src2 & 0xFF);
    u32Min1 = (u32Src1 & 0xFF00) < (u32Src2 & 0xFF00)? (u32Src1 & 0xFF00):(u32Src2 & 0xFF00);
    u32Min2 = (u32Src1 & 0xFF0000) < (u32Src2 & 0xFF0000)? (u32Src1 & 0xFF0000):(u32Src2 & 0xFF0000);
    u32Min3 = (u32Src1 & 0xFF000000) < (u32Src2 & 0xFF000000)? (u32Src1 & 0xFF000000):(u32Src2 & 0xFF000000);

    return u32Min0|u32Min1|u32Min2|u32Min3;
}
static __inline  u32 _add4(u32 u32Src1, u32 u32Src2)
{
    u32 u32ret=0;

    u32ret |= ((BYTE0(u32Src1) + BYTE0(u32Src2)) & 0xFF);
    u32ret |= ((BYTE1(u32Src1) + BYTE1(u32Src2)) & 0xFF) << 8;
    u32ret |= ((BYTE2(u32Src1) + BYTE2(u32Src2)) & 0xFF) << 16;
    u32ret |= ((BYTE3(u32Src1) + BYTE3(u32Src2)) & 0xFF) << 24;
    
    return u32ret;
}
static __inline  u32 _saddu4(u32 u32Src1, u32 u32Src2)
{
    u32 u32sum0;
    u32 u32sum1;
    u32 u32sum2;
    u32 u32sum3;

    u32sum0 = (BYTE0(u32Src1) + BYTE0(u32Src2));
    u32sum0 = u32sum0 > 255? 255: u32sum0;

    u32sum1 = (BYTE1(u32Src1) + BYTE1(u32Src2));
    u32sum1 = u32sum1 > 255? 255: u32sum1;

    u32sum2 = (BYTE2(u32Src1) + BYTE2(u32Src2));
    u32sum2 = u32sum2 > 255? 255: u32sum2;
    
    u32sum3 = (BYTE3(u32Src1) + BYTE3(u32Src2));
    u32sum3 = u32sum3 > 255? 255: u32sum3;

    return (u32sum3 << 24)|(u32sum2 << 16)|(u32sum1 << 8)| u32sum0;
}
static __inline  u32 _sub4(u32 u32Src1, u32 u32Src2)
{
    u32 u32ret=0;

    u32ret |= ((BYTE0(u32Src1) - BYTE0(u32Src2)) & 0xFF);
    u32ret |= ((BYTE1(u32Src1) - BYTE1(u32Src2)) & 0xFF) << 8;
    u32ret |= ((BYTE2(u32Src1) - BYTE2(u32Src2)) & 0xFF) << 16;
    u32ret |= ((BYTE3(u32Src1) - BYTE3(u32Src2)) & 0xFF) << 24;
    
    return u32ret;
}
static __inline  u32 _avgu4(u32 u32Src1, u32 u32Src2)
{
    u32 u32ret=0;
    
    u32ret |= (((BYTE0(u32Src1) + BYTE0(u32Src2) + 1) >> 1) & 0xFF);
    u32ret |= (((BYTE1(u32Src1) + BYTE1(u32Src2) + 1) >> 1) & 0xFF) << 8;
    u32ret |= (((BYTE2(u32Src1) + BYTE2(u32Src2) + 1) >> 1) & 0xFF) << 16;
    u32ret |= (((BYTE3(u32Src1) + BYTE3(u32Src2) + 1) >> 1) & 0xFF) << 24;
    
    return u32ret;
}
static __inline  s64 _dpack2(u32 u32Src1, u32 u32Src2)
{
    u32 u32A3, u32A2;
    u32A3 = _packh2(u32Src1, u32Src2);
    u32A2 = _pack2(u32Src1, u32Src2);
    return _itoll(u32A3, u32A2);
}
static __inline  u32 _abs(l32 l32Src1)
{
    return ABS(l32Src1);
}
static __inline u32 _bitc4(u32 u32Src1)
{
    u32 u32BitCnt;
    l32 l32Index;
    l32 l32Unit;
    
    u32BitCnt = 0;
    
    l32Unit = 0x01000000;
    for(l32Index = 31; l32Index >= 0; l32Index--)
    {
        if(u32Src1 & (1 << l32Index))
        {
            u32BitCnt += l32Unit;
        }

        if((l32Index % 8) == 0)
        {
            l32Unit >>= 8;
        }
    }

    return u32BitCnt;
}
static __inline u32 _sshvr(l32 l32Src1, l32 l32a)
{
    if(l32a > 31)
    {
        l32a=31;
    }
    if(l32a < -31)
    {
        l32a=-31;
    }

    if(l32a>=0)
    {
        return (l32Src1 >> l32a);
    }
    else
    {
        //饱和左移
        l32a = ABS(l32a);
        
        if( (l32Src1 > 0) && (l32Src1 >> (32 - l32a)) )
        {
            return 0x7FFFFFFF;
        }
        if( (l32Src1 < 0) && ((l32Src1 >> (32 - l32a)) != -1) )
        {
            return 0x80000000;
        }
        
        return (l32Src1 << l32a);
    }
}
static __inline u32 _sshvl(l32 l32Src1, l32 l32a)
{
    if(l32a > 31)
    {
        l32a = 31;
    }
    if(l32a < -31)
    {
        l32a = -31;
    }
    
    if(l32a < 0)
    {
        return (l32Src1 >> ABS(l32a));
    }
    else
    {
        //饱和左移
        
        if( (l32Src1 > 0) && (l32Src1 >> (32 - l32a)) )
        {
            return 0x7FFFFFFF;
        }
        if( (l32Src1 < 0) && ((l32Src1 >> (32 - l32a)) != -1) )
        {
            return 0x80000000;
        }
        
        return (l32Src1 << l32a);
    }
}
static __inline u32 _shlmb(u32 u32Src1, u32 u32Src2)
{
    return (u32Src2 << 8) | (u32Src1 >> 24);
}    
static __inline u32 _shrmb(u32 u32Src1, u32 u32Src2)
{
    return (u32Src2 >> 8) | ((u32Src1 & 0xFF) << 24);
}    
//下面的函数可能造成编译警告
//      warning C4142: benign redefinition of type
//      warning C4028: formal parameter 1 different from declaration
//      warning C4028: formal parameter 2 different from declaration
//
// 这是因为标准C运行时库提供了循环左移的intrinsic实现,
// 但是我们需要使用自己的实现确保模拟TI DSP的对应指令.
// 因此请忽略这样的警告.
/*
static __inline u32 _rotl(u32 u32Src1, u32 u32a)
{
    u32 u32Ret;
    u32Ret =  u32Src1 << u32a;
    u32Ret |= u32Src1 >> (32 - u32a);
    return u32Ret;
}
*/
static __inline u32 _swap4(u32 u32Src1)
{
    u32 u32Ret=0;
    u32Ret |= (u32Src1 << 8) & 0xFF00FF00;
    u32Ret |= (u32Src1 >> 8) & 0x00FF00FF;
    return u32Ret;
}


//下面的内存访问宏必须兼容出现在等号左边和右边的情况
//	使用内联函数返回C++中的引用类型
//	或者使用宏定义配合C中的指针访问都可以实现
#define _mem4(p) (*(u32 *)(p))
#define _amem4(p) (*(u32 *)(p))
#define _amem4_const(p) (*(u32 *)(p))
#define _mem4_const(p) (*(u32 *)(p))

#define _mem8(p) (*(s64 *)(p))
#define _amem8(p) (*(s64 *)(p))
#define _memd8(p) (*(x64 *)(p))
#define _amemd8(p) (*(x64 *)(p))
#define _amemd8_const(p) (*(x64 *)(p))
#define _memd8_const(p) (*(x64 *)(p))
#define _amem8_const(p) (*(s64 *)(p))
#define _mem8_const(p) (*(s64 *)(p))


//由于输入参数d可能是double(也是64bit)类型,因此位操作前转换为x64类型
#define _lo(d) (u32)(  (*(x64*)(&d))       & 0xFFFFFFFFL)
#define _hi(d) (u32)(( (*(x64*)(&d)) >>32) & 0xFFFFFFFFL)

//VC环境下的右移操作在移动量大于等于32时未定义(由于是32bit移位硬件)
//因此使用下面的inline函数模拟TI DSP右移32位时结果为0的操作
#define _SHRU(a,b)  _shru_ti(a,b)

typedef unsigned __int64 u40;

#define _loll(d) (u32)(  (*(u64*)(&d))       & 0xFFFFFFFFL)
#define _hill(d) (u32)(( (*(u64*)(&d)) >>32) & 0xFFFFFFFFL)
#define _mem2(p) (*(u16 *)(p))
#define _amem2(p) (*(u16 *)(p))

/*=========================================================================
                                   新增
==========================================================================*/

#define U64_BYTE16(u64A) ((u64)(u64A) & 0xFFFF)
#define U64_BYTE32(u64A) ((u64)(u64A>>16) & 0xFFFF)
#define U64_BYTE48(u64A) ((u64)(u64A>>32) & 0xFFFF)
#define U64_BYTE64(u64A) ((u64)(u64A>>48) & 0xFFFF)

static __inline u64 _dadd2(u64 u64Src1, u64 u64Src2)
{
    u64 u64Res;

    u64Res = ((U64_BYTE64(u64Src1) + U64_BYTE64(u64Src2)) & 0xFFFF) << 16;    
    u64Res = (u64Res | (((U64_BYTE48(u64Src1) + U64_BYTE48(u64Src2)) & 0xFFFF))) << 16;
    u64Res = (u64Res | ((U64_BYTE32(u64Src1) + U64_BYTE32(u64Src2)) & 0xFFFF)) << 16;
    u64Res = (u64Res | ((U64_BYTE16(u64Src1) + U64_BYTE16(u64Src2)) & 0xFFFF));
    
    return u64Res;
}

static __inline u64 _dsub2(u64 u64Src1, u64 u64Src2)
{
    u64 u64Res;

    u64Res = ((U64_BYTE64(u64Src1) - U64_BYTE64(u64Src2)) & 0xFFFF) << 16;    
    u64Res = (u64Res | (((U64_BYTE48(u64Src1) - U64_BYTE48(u64Src2)) & 0xFFFF))) << 16;
    u64Res = (u64Res | ((U64_BYTE32(u64Src1) - U64_BYTE32(u64Src2)) & 0xFFFF)) << 16;
    u64Res = (u64Res | ((U64_BYTE16(u64Src1) - U64_BYTE16(u64Src2)) & 0xFFFF));

    return u64Res;
}

static __inline  u64 _unpkbu4(u32 u32Src1)
{
    u64 u64Ret = 0;

    u64Ret = (u64Ret | _extu(u32Src1, 0, 24)) << 16;
    u64Ret = (u64Ret | _extu(u32Src1, 8, 24)) << 16;
    u64Ret = (u64Ret | _extu(u32Src1, 16, 24)) << 16;
    u64Ret = (u64Ret | _extu(u32Src1, 24, 24));

    return u64Ret;
}

static __inline u64 _mpyu4ll(u32 u32Src1, u32 u32Src2)
{
    u64 u64ret=0;

    u64ret |= ((BYTE3(u32Src1) * BYTE3(u32Src2)) & 0xFFFF);
    u64ret <<= 16;
    u64ret |= ((BYTE2(u32Src1) * BYTE2(u32Src2)) & 0xFFFF);
    u64ret <<= 16;
    u64ret |= ((BYTE1(u32Src1) * BYTE1(u32Src2)) & 0xFFFF);
    u64ret <<= 16;
    u64ret |= ((BYTE0(u32Src1) * BYTE0(u32Src2)) & 0xFFFF);

    return u64ret;
}

#define U64_BYTE_32(u64A) ((u64)(u64A) & 0xFFFFFFFF)
#define U64_BYTE_64(u64A) ((u64)(u64A>>32) & 0xFFFFFFFF)
static __inline u64 _dadd(u64 u64Src1, u64 u64Src2)
{
    u64 u64Res;

    u64Res = ((U64_BYTE_64(u64Src1) + U64_BYTE_64(u64Src2)) & 0xFFFFFFFF) << 32;    
    u64Res = u64Res | ((U64_BYTE_32(u64Src1) + U64_BYTE_32(u64Src2)) & 0xFFFFFFFF);

    return u64Res;
}

static __inline u64 _unpkhu2(u32 u32Src1)
{
    u64 u64Ret = 0;
    u64Ret = (u64Ret | _extu(u32Src1, 0, 16)) << 32;
    u64Ret = (u64Ret | _extu(u32Src1, 16, 16));
    return u64Ret;
}

static __inline float _rsqrsp (float src)
{
    float f32Res;
    f32Res = sqrt(src);
    return 1.0f / f32Res;
}

static __inline float _rcpsp (float src)
{
    return 1.0f / src;
}
static __inline u64 _dshru(u64 u64Src1, u32 u32Src2)
{
    u32 u32H = (u64Src1 >> 32) & 0xFFFFFFFFL;
    u32 u32L = (u64Src1 & 0xFFFFFFFFL);
    
    u32H >>= u32Src2;
    u32L >>= u32Src2;

    u64Src1 = u32H;
    return (u64Src1 << 32) | ((u64)u32L);
}

static __inline u64 _dshr2(u64 u64Src1, u32 u32Src2)
{
    u32 u32H = (u64Src1 >> 32) & 0xFFFFFFFFL;
    u32 u32L = (u64Src1 & 0xFFFFFFFFL);
    
    u32H = _shr2(u32H, u32Src2);
    u32L = _shr2(u32L, u32Src2);

    u64Src1 = u32H;
    return (u64Src1 << 32) | ((u64)u32L);
}

static __inline u64 _dshl(u64 u64Src1, u32 u32Src2)
{
    u32 u32H = (u64Src1 >> 32) & 0xFFFFFFFFL;
    u32 u32L = (u64Src1 & 0xFFFFFFFFL);
    
    u32H <<= u32Src2;
    u32L <<= u32Src2;

    u64Src1 = u32H;
    return (u64Src1 << 32) | ((u64)u32L);
}

static __inline long long _mpy2ll (int src1, int src2)
{
    u64 u64ret = 0;

    u64ret |= ((_ext(src1, 0, 16) * _ext(src2, 0, 16)) & 0xFFFFFFFF);
    u64ret <<= 32;
    u64ret |= ((_ext(src1, 16, 16) * _ext(src2, 16, 16)) & 0xFFFFFFFF);

    return u64ret;
}

static __inline u64 _mpyu2(u32 u32Src1,u32 u32Src2)
{
    u32 u32L = (u32Src1 & 0xFFFF) * (u32Src2 & 0xFFFF);
    u32 u32H = (u32Src1 >> 16) * (u32Src2 >> 16);
    u64 u64Ret;
    u64Ret = u32H;
    return (u64Ret << 32) | ((u64)u32L);
}

static __inline u64 _dpackhl2(u64 u64Src1, u64 u64Src2)
{
    u32 u32H1 = (u64Src1 >> 32) & 0xFFFFFFFFL;
    u32 u32L1 = (u64Src1 & 0xFFFFFFFFL);

    u32 u32H2 = (u64Src2 >> 32) & 0xFFFFFFFFL;
    u32 u32L2 = (u64Src2 & 0xFFFFFFFFL);

    u32H2 = _packhl2(u32H1, u32H2);
    u32L2 = _packhl2(u32L1, u32L2);

    u64Src1 = u32H2;
    return (u64Src1 << 32) | ((u64)u32L2);
}

static __inline u64 _dpackh2(u64 u64Src1, u64 u64Src2)
{
    u32 u32H1 = (u64Src1 >> 32) & 0xFFFFFFFFL;
    u32 u32L1 = (u64Src1 & 0xFFFFFFFFL);

    u32 u32H2 = (u64Src2 >> 32) & 0xFFFFFFFFL;
    u32 u32L2 = (u64Src2 & 0xFFFFFFFFL);

    u32H2 = _packh2(u32H1, u32H2);
    u32L2 = _packh2(u32L1, u32L2);

    u64Src1 = u32H2;
    return (u64Src1 << 32) | ((u64)u32L2);
}

static __inline u64 _dpackl2(u64 u64Src1, u64 u64Src2)
{
    u32 u32H1 = (u64Src1 >> 32) & 0xFFFFFFFFL;
    u32 u32L1 = (u64Src1 & 0xFFFFFFFFL);

    u32 u32H2 = (u64Src2 >> 32) & 0xFFFFFFFFL;
    u32 u32L2 = (u64Src2 & 0xFFFFFFFFL);

    u32H2 = _pack2(u32H1, u32H2);
    u32L2 = _pack2(u32L1, u32L2);

    u64Src1 = u32H2;
    return (u64Src1 << 32) | ((u64)u32L2);
}
static __inline u32 _dcmpgtu4 (s64 src1, s64 src2)
{
    u32 u32lo = _cmpgtu4(_loll(src1), _loll(src2));
    u32 u32hi = _cmpgtu4(_hill(src1), _hill(src2));
    return (u32hi << 4)|(u32lo);
}
static __inline s64 _dxpnd4 (u32 u32src)
{
    return _itoll(_xpnd4(u32src >> 4), _xpnd4(u32src));
}


/*=========================================================================
                                   新增
==========================================================================*/

#endif

#endif //#ifndef _INTRINSIC_H_


