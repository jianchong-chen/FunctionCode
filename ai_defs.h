/******************************************************************************
模块名      ： ai_defs
******************************************************************************/
#ifndef __AIDEFS__H__
#define __AIDEFS__H__

#ifdef __cplusplus
extern "C" {
#endif

//==========================类型定义==========================
#ifdef NEED_TYPE_BOOL32
    typedef int BOOL32;
#endif

typedef char            s8;            
typedef short           s16;
typedef unsigned char   u8;   
typedef unsigned short  u16;	
typedef float           f32;
typedef double          d64;

#if defined _TMS320C6400 || defined(ARCH_C667X) || defined (ARCH_C674X) //TI DM64X平台宏定义
    typedef double             x64;
    typedef long long		   s64;
    typedef unsigned long long u64;
    typedef int                l32;
    typedef unsigned int       u32;
	
	#define INLINE inline
	#define FAST_CALL 
	#define EXPORT 
    #define RESTRICT restrict
#elif defined(ARCH_X86_LINUX) || defined(ARCH_ATOM_LINUX) || defined(ARCH_ARM_LINUX)  || defined(ARCH_POWERPC_LINUX)//(x86及ARM平台)Linux系统
    typedef unsigned long long x64;
    typedef long long		   s64;
    typedef unsigned long long u64;
    typedef long               l32;
    typedef unsigned long      u32;
	
	#define INLINE inline
	#define FAST_CALL __attribute__((fastcall))
	#define EXPORT 
    #define RESTRICT 
#else //目前仅限Windows系统
    typedef unsigned __int64   x64;
    typedef __int64			   s64;
    typedef unsigned __int64   u64;
    typedef long               l32;
    typedef unsigned long      u32;
	
	#define INLINE __inline
	#define FAST_CALL __fastcall
	//Dll导出支持
    #define EXPORT __declspec(dllexport)
    #define RESTRICT 
#endif

//Ubuntu系统定义
#ifdef ARCH_ARM_UBUNTU
    typedef signed char     s8;         
    typedef signed short    s16;
	  typedef signed __int64	s64;
#endif

// s8在有的平台（如c665x）需要显式定义为有符号的
#ifdef NEED_TYPE_SIGNED_CHAR
    typedef signed char     s8;
#endif

#ifndef TRUE
    #define TRUE    1
#endif

#ifndef FALSE
    #define FALSE   0
#endif

#ifndef NULL
  #ifdef  __cplusplus
    #define NULL    0
  #else
    #define NULL    ((void *)0)
  #endif
#endif

#ifndef ARCH_ARM_WINCE
//typedef int BOOL;
#if !defined(TYPE_BOOL) && !defined(__INCvxTypesh)
    typedef int BOOL, *PBOOL;
    #define TRUE 1
    #define FALSE 0
    #define TYPE_BOOL
#endif /* BOOL */

#endif

//==========================公共宏定义==========================
#ifndef MIN
#define MIN(a, b)   ((a) > (b) ? (b) : (a))
#endif
#ifndef MAX
#define MAX(a, b)   ((a) > (b) ? (a) : (b))
#endif
#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : (-(x)))
#endif

//==========================公共结构体定义======================
#define MAX_NUM_LINE_POINT  64        //折线的点的最多个数
#define MAX_NUM_POLYGON_POINT  64     //多边形的点的最多个数
typedef struct
{
    l32 l32X;
    l32 l32Y;
}TPoint;

typedef struct
{
    f32 f32X;
    f32 f32Y;
}TFPoint;

typedef struct                      //浮点折线
{
    l32 l32PointNum;
    TPoint atPoint[MAX_NUM_LINE_POINT];
}TLine;

typedef struct                      //浮点折线
{
    l32 l32FPointNum;
    TFPoint atFPoint[MAX_NUM_LINE_POINT];
}TFLine;

//不包含右边界，及下边界
//l32Left <= x < l32Top + l32Width;
//l32Top <= x < l32Top + l32Height;
typedef struct
{
    l32 l32Left;  //x坐标
    l32 l32Top;   //y坐标
    l32 l32Width;
    l32 l32Height;
}TRect;

typedef struct 
{
    l32 l32PointNum;    
    TPoint atPoint[MAX_NUM_POLYGON_POINT];
}TPolygon;

typedef struct 
{
    l32 l32FPointNum;    
    TFPoint atFPoint[MAX_NUM_POLYGON_POINT];
}TFPolygon;

typedef struct
{
    l32 l32Width;   //宽度
    l32 l32Height;  //高度
}TSize;

//定义同图像属性，类型，相关的常量
enum EImageAttribute
{
    //图像属性：Attribute
    //FRAME_I_FORMAT = INTERLACE_CAPTURE
    //FRAME_FORMAT = 0
    //FIELD_FORMAT = INTERLACE_CAPTURE | FIELD_STORE
    FIELD_STORE = (1 << 0),             //0bit:隔行存放
    INTERLACE_CAPTURE = (1 << 1),       //1bit:隔行采集    
	FLIP_V = (1 << 2),                  //2bit:垂直翻转
	FLIP_H = (1 << 3)                   //3bit:水平翻转
};

#define ATTRIBUTE_FRAME_I_FORMAT (INTERLACE_CAPTURE)
#define ATTRIBUTE_FRAME_FORMAT (0)
#define ATTRIBUTE_FIELD_FORMAT (INTERLACE_CAPTURE | FIELD_STORE)

#define MAKE_FOURCC(a, b, c, d) (((u32)(a) << 24) | ((u32)(b) << 16) | ((u32)(c) << 8) | (u32)(d))

//可以图像类型为，Attribute + Type
typedef enum 
{ 
    //无法用FOURCC表示的类型定义在此处
    AI_RGB16,
    AI_RGB24,
    AI_RGB = AI_RGB24,
    AI_RGB32,
    AI_HSV,
    AI_Y,

    //=======================================================
    //YUV数据格式，同存储方式紧密结合
    //-------------------------------
    //YUV444: Y Plane: U Plane: V Plane
    //      Horizontal      Vertical
    //Y         1               1
    //UV        1               1
    AI_Y444 = MAKE_FOURCC('Y', '4', '4', '4'),
    //-------------------------------
    //UYVY(Y422): YUV 4:2:2 
    //      [U0 Y0 V0 Y1][U2 Y2 V2 Y3]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               1
    AI_UYVY = MAKE_FOURCC('U', 'Y', 'V', 'Y'),
    AI_Y422 = MAKE_FOURCC('Y', '4', '2', '2'),
    //-------------------------------
    //YUYV: YUV 4:2:2 
    //      [Y0 U0 Y1 V0][Y2 U2 Y3 V2]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               1
    AI_YUYV = MAKE_FOURCC('Y', 'U', 'Y', 'V'),
    //-------------------------------
    //YVYU: YUV 4:2:2 
    //      [Y0 V0 Y1 U0][Y2 V2 Y3 U2]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               1
    AI_YVYU = MAKE_FOURCC('Y', 'V', 'Y', 'U'),
    //-------------------------------
    //YV16: YUV422:
    //      Y[W,H] U[W/2, H] V[W/2, H]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               1
    AI_YV16 = MAKE_FOURCC('Y', 'V', '1', '6'),
    //-------------------------------
    //I422: YUV422:
    //      Y[N,M] U[N/2,M] V[N/2, M]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               1
    AI_I422 = MAKE_FOURCC('I', '4', '2', '2'),
    //-------------------------------
    //I420: YUV420:
    //      Y[N,M] U[N/2,M/2] V[N/2, M/2]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               2
    AI_I420 = MAKE_FOURCC('I', '4', '2', '0'),
    //-------------------------------
    //YV12: YUV420:
    //      Y[N,M] V[N/2,M/2] U[N/2, M/2]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               2
    AI_YV12 = MAKE_FOURCC('Y', 'V', '1', '2'),
    //-------------------------------
    //NV12: YUV420:
    //      Y[N,M] UV[N/2,M/2]
    //      [U0,V0] [U1,V1]
    //      Horizontal      Vertical
    //Y         1               1
    //UV        2               2
    AI_NV12 = MAKE_FOURCC('N', 'V', '1', '2')
}EImageType;

typedef  struct  
{
    l32 l32Width;
    l32 l32Height;
    l32 l32Stride;
    void *pvBuffer;
}TPlane;

typedef struct
{
	u32 u32Attribute;      //图像属性,具体内容参见EImageAttribute
	u32 u32Type;           //图像类型,具体内容参见EImageType
    //如果为YUV交织模式:
    //      atPlane[0]:存放交织数据
    //      atPlane[1,2,3]为0,无效
    //如果为YUV Plane模式，
    //      atPlane[0]:存放Y数据
    //      atPlane[1]:存放U数据
    //      atPlane[2]:存放V数据
    //      atPlane[3]为0,无效
    //如果为Y Plane UV交织存放模式：
    //      atPlane[0]:存放Y数据
    //      atPlane[1]:存放交织的UV数据
    //      atPlane[2,3]为0,无效
    //如果为RGB+阿尔法通道则：
    //      atPlane[0]:存放R数据
    //      atPlane[1]:存放G数据
    //      atPlane[2]:存放B数据
    //      atPlane[3]:存放阿尔法通道数据
	TPlane atPlane[4];
	u32 u32Reserved;
}TImage;

//增加时间戳的TImage
typedef struct
{
    s64 s64TimeStamp; //绝对时间ms   
    TImage tImg;	
}TFrame;

//提供宏定义设置TImage图像格式:
//定义RGB24连续buffer图像结构体
#define SET_TIMAGE_RGB24_CONTINUOUS(pImageIn, width, height, Buffer) \
UtiSetTImage(pImageIn, 0, AI_RGB24, width, height, Buffer, NULL, NULL, NULL, 0, 0, 0, 0);
//定义RGB32连续buffer图像结构体
#define SET_TIMAGE_RGB32_CONTINUOUS(pImageIn, width, height, Buffer) \
UtiSetTImage(pImageIn, 0, AI_RGB32, width, height, Buffer, NULL, NULL, NULL, 0, 0, 0, 0);
//定义I420连续buffer图像结构体
#define SET_TIMAGE_I420_CONTINUOUS(pImageIn, width, height, Buffer) \
UtiSetTImage(pImageIn, 0, AI_I420, width, height, Buffer, NULL, NULL, NULL, 0, 0, 0, 0);
//定义YV12连续buffer图像结构体
#define SET_TIMAGE_YV12_CONTINUOUS(pImageIn, width, height, Buffer) \
UtiSetTImage(pImageIn, 0, AI_YV12, width, height, Buffer, NULL, NULL, NULL, 0, 0, 0, 0);
//定义NV12连续buffer图像结构体
#define SET_TIMAGE_NV12_CONTINUOUS(pImageIn, width, height, Buffer) \
UtiSetTImage(pImageIn, 0, AI_NV12, width, height, Buffer, NULL, NULL, NULL, 0, 0, 0, 0);
//定义YV16连续buffer图像结构体
#define SET_TIMAGE_YV16_CONTINUOUS(pImageIn, width, height, Buffer) \
UtiSetTImage(pImageIn, 0, AI_YV16, width, height, Buffer, NULL, NULL, NULL, 0, 0, 0, 0);

//定义I420不连续buffer图像结构体,buffer Stride等同于各个分量宽度
#define SET_TIMAGE_I420_SEPARATED(pImageIn, width, height, Buffer0, Buffer1, Buffer2) \
UtiSetTImage(pImageIn, 0, AI_I420, width, height, Buffer0, Buffer1, Buffer2, NULL, 0, 0, 0, 0);
//定义YV12不连续buffer图像结构体,buffer Stride等同于各个分量宽度
#define SET_TIMAGE_YV12_SEPARATED(pImageIn, width, height, Buffer0, Buffer1, Buffer2) \
UtiSetTImage(pImageIn, 0, AI_YV12, width, height, Buffer0, Buffer1, Buffer2, NULL, 0, 0, 0, 0);
//定义NV12不连续buffer图像结构体,buffer Stride等同于各个分量宽度
#define SET_TIMAGE_NV12_SEPARATED(pImageIn, width, height, Buffer0, Buffer1, Buffer2) \
UtiSetTImage(pImageIn, 0, AI_NV12, width, height, Buffer0, Buffer1, Buffer2, NULL, 0, 0, 0, 0);
//定义YV16不连续buffer图像结构体,buffer Stride等同于各个分量宽度
#define SET_TIMAGE_YV16_SEPARATED(pImageIn, width, height, Buffer0, Buffer1, Buffer2) \
UtiSetTImage(pImageIn, 0, AI_YV16, width, height, Buffer0, Buffer1, Buffer2, NULL, 0, 0, 0, 0);

//定义I420不连续buffer图像结构体,buffer Stride自己设定
#define SET_TIMAGE_I420_SEPARATED_STEP(pImageIn, width, height, Buffer0, Buffer1, Buffer2, Stride0, Stride1, Stride2) \
UtiSetTImage(pImageIn, 0, AI_I420, width, height, Buffer0, Buffer1, Buffer2, NULL, Stride0, Stride1, Stride2, 0);
//定义YV12不连续buffer图像结构体,buffer Stride自己设定
#define SET_TIMAGE_YV12_SEPARATED_STEP(pImageIn, width, height, Buffer0, Buffer1, Buffer2, Stride0, Stride1, Stride2) \
UtiSetTImage(pImageIn, 0, AI_YV12, width, height, Buffer0, Buffer1, Buffer2, NULL, Stride0, Stride1, Stride2, 0);
//定义NV12不连续buffer图像结构体,buffer Stride自己设定
#define SET_TIMAGE_NV12_SEPARATED_STEP(pImageIn, width, height, Buffer0, Buffer1, Buffer2, Stride0, Stride1, Stride2) \
UtiSetTImage(pImageIn, 0, AI_NV12, width, height, Buffer0, Buffer1, Buffer2, NULL, Stride0, Stride1, Stride2, 0);
//定义YV16不连续buffer图像结构体,buffer Stride自己设定
#define SET_TIMAGE_YV16_SEPARATED_STEP(pImageIn, width, height, Buffer0, Buffer1, Buffer2, Stride0, Stride1, Stride2) \
UtiSetTImage(pImageIn, 0, AI_YV16, width, height, Buffer0, Buffer1, Buffer2, NULL, Stride0, Stride1, Stride2, 0);

/*====================================================================
函数名      ：  UtiSetTImage
功能        ：	设置图像结构体
引用全局变量：	无
输入参数说明：  Out:   *pImageIn:       图像结构体
                In:     u32Attribute    图像属性
                In:     u32Type         图像类型
                In:     u32Width        图像宽度
                In:     u32Height       图像高度
                In:     pvBuffer0       buffer0地址
                In:     pvBuffer1       buffer1地址
                In:     pvBuffer2       buffer2地址
                In:     pvBuffer3       buffer3地址
                In:     u32Stride0      buffer0步长
                In:     u32Stride1      buffer1步长
                In:     u32Stride2      buffer2步长
                In:     u32Stride3      buffer3步长
For example:
                1:连续YUV，stride等于宽度
                UtiSetTImage(pImageSrc, 0, I420, 352, 288, 
                          u8YUVBuf, 0, 0, 0,
                          0, 0, 0, 0);
                2:VU buffer连续，Y buffer 独立
                UtiSetTImage(pImageSrc, 0, YV16,  352, 288, 
                          u8YBuf, u8UBuf, 0, 0,
                          720, 352, 0, 0);
                3:YU buffer连续， V buffer独立 场格式,水平翻转
                UtiSetTImage(pImageSrc, INTERLACE | FIELD | FLIP_H, I420,  352, 288, 
                          u8YBuf, 0, u8VBuf, 0,
                          1024, 0, 720, 0);
回值说明  ：	无
=============================================================================*/
EXPORT void UtiSetTImage(TImage *pImageIn, 
                         u32 u32Attribute, u32 u32Type, 
                         l32 l32Width, l32 l32Height,
                         void *pvBuffer0, void *pvBuffer1, void *pvBuffer2,void *pvBuffer3,
                         l32 l32Stride0, l32 l32Stride1, l32 l32Stride2, l32 l32Stride3);
                      
/*====================================================================
函数名      ：  UtiIsTImageLegal
功能        ：  判断图像结构体是否合法
引用全局变量：  无
输入参数说明：  in:   pImageIn:       图像结构体
                  in:   u32WDivisor0, u32HDivisor0, Plane0宽高约数
                  in:   u32WDivisor1, u32HDivisor1, Plane1,Plane2,Plane3宽高约数
				  in:   u32Align0:      Plane0对齐参数
                  in:   u32Align1:      Plane1对齐参数
                  in:   u32Align2:      Plane2对齐参数
                  in:   u32Align3:      Plane3对齐参数
回值说明  ：	 u32Result：
                 32                 --->               0
                 .... .... .... ..xx xxxx xxxx xxxx xxxx
Field[3-0]: 表示Plane3是否正常：
            bit0：buffer为空
            bit1：buffer没有满足对齐要求
            bit2：buffer宽高超出最大范围或没有满足约数要求
            bit3：buffer步长小于宽度
Field[7-4]: 表示Plane2是否正常：
            bit4：buffer为空
            bit5：buffer没有满足对齐要求
            bit6：buffer宽高超出最大范围或没有满足约数要求
            bit7：buffer步长小于宽度
Field[11-8]: 表示Plane1是否正常：
            bit8：buffer为空
            bit9：buffer没有满足对齐要求
            bit10：buffer宽高超出最大范围或没有满足约数要求
            bit11：buffer步长小于宽度
Field[15-12]: 表示Plane0是否正常：
            bit12：buffer为空
            bit13：buffer没有满足对齐要求
            bit14：buffer宽高超出最大范围或没有满足约数要求
            bit15：buffer步长小于宽度
Field[17-16]:
            bit16：图像属性组合错误
            bit17：图像类型错误
=============================================================================*/                      
EXPORT u32 UtiIsTImageLegal(const TImage *pImageIn, 
                            u32 u32WDivisor0, u32 u32HDivisor0, 
                            u32 u32WDivisor1, u32 u32HDivisor1, 
                            u32 u32Align0, u32 u32Align1, u32 u32Align2, u32 u32Align3);

//==========================错误码定义==========================
typedef enum
{
    EStatus_Success = 0,                             // 成功
    EStatus_InvalidParameter = 0x80000001,           // 参数错误
    EStatus_OutOfMemory = 0x80000002,                // 内存分配失败
    EStatus_InsufficientBuffer = 0x80000003,         // buffer大小不够
    EStatus_GenericError = 0x80000004,               // 特定错误
    EStatus_Undefined = 0x80000005                   // 未知错误
}EStatus;
     
//==========================告警类型============================         
//与用户设置方案一致
typedef enum
{

    //普通周界告警功能
    ALARM_CROSS_LINE = 4,           //过线告警
    ALARM_CROSS_WALL,           //翻墙告警
    ALARM_CROSS_AREA,           //越区域告警
    ALARM_WANDER,               //徘徊告警

    //行为分析告警
    ALARM_RABBLE = 4000,               //聚众告警
    ALARM_RUNNING,              //奔跑告警
    ALARM_FIGHTING,              //打架告警
    ALARM_AREAGUARD_IN,          //区域看防进入告警
    ALARM_AREAGUARD_OVERTIME,    //区域看防逗留超时告警
    ALARM_AREAGUARD_OUT,         //区域看防出去告警
    ALARM_RISE,                  //起身检测告警
    ALARM_ONDUTYDETECT_IN,          //值岗检测进入告警
    ALARM_ONDUTYDETECT_OVERTIME,    //值岗检测逗留超时告警
    ALARM_ONDUTYDETECT_OUT,         //值岗检测出去告警
	ALARM_SPEECHDETECT,             //声音检测异常告警
    
    //视质轮询告警
	ALARM_BLUR = 5000,                 //模糊
	ALARM_BRIGHTHIGH,            //过亮
	ALARM_BRIGHTLOW,             //过暗
	ALARM_BRIGHTCON,             //对比度异常
	ALARM_FREEZE,                //画面冻结
	ALARM_HUE,                   //色偏
	ALARM_INTERFERENCE,          //干扰
	ALARM_PTZSTATIC,             //PTZ检测静止
	ALARM_PTZLEFT,               //PTZ检测向左
	ALARM_PTZRIGHT,              //PTZ检测向右
	ALARM_PTZUP,                 //PTZ检测向上
	ALARM_PTZDOWN,               //PTZ检测向下
	ALARM_PTZUPLEFT,             //PTZ检测向左上
	ALARM_PTZUPRIGHT,            //PTZ检测向右上
	ALARM_PTZDOWNLEFT,           //PTZ检测向左下
	ALARM_PTZDOWNRIGHT,          //PTZ检测向右下
	ALARM_PTZZOOMIN,             //PTZ检测ZoomIn
	ALARM_PTZZOOMOUT,            //PTZ检测ZoomOut
	ALARM_SHAKE,                 //抖动
	ALARM_SHELTER,               //遮挡
	ALARM_SIGLOSS,               //信号丢失
	ALARM_SCENECHANGE            //场景切换
}EAlarmType;

#ifdef __cplusplus
}
#endif

#endif
