#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "opt_inc.h"

#include "ext_blob.h"
#include <string.h>


/*=============================================================================
函数名    ：CAntImage::CAntImage()
功能      ：类CAntImage的构造函数
算法实现  ：无
参数说明  ：vw  
vh  
返回值说明：无
其他说明  ：无
=============================================================================*/
CAntImage::CAntImage(int vw, int vh):m_width(vw),m_height(vh),m_nChannels(3)
{
    int i;
    //allocates memory for accompanying images
    m_pTempf1 = new int[vw*vh];

    m_pTempMask = new unsigned char[vw*vh];
    memset(m_pTempMask, 0, vw * vh*sizeof(unsigned char));

    for(i=0;i<MAX_NUM_BLOBS;i++)
    {
        m_blobs[i] = new AntBlob();
    }

    m_nBlobs = 0;
    m_minTargetSizeX = 3; //4 by 4 pixels
    m_minTargetSizeY = 3; //4 by 4 pixels
    m_x_dilate_size = 8; //size of x direction dilation
    m_y_dilate_size = 10; //size of y direction dilation
}

/*=============================================================================
函数名    ：CAntImage::~CAntImage()
功能      ：类CAntImage的析构函数
算法实现  ：无
参数说明  ：无
返回值说明：无
其他说明  ：无
=============================================================================*/
CAntImage::~CAntImage()
{
    FreeAllMemory();
}

/*=============================================================================
函数名    ：CAntImage::FreeAllMemory()
功能      ：释放所有的内存
算法实现  ：无
参数说明  ：无
返回值说明：无
其他说明  ：无
=============================================================================*/
void CAntImage::FreeAllMemory()
{
    //allocates memory for accompanying images
    delete [] m_pTempf1;

    delete [] m_pTempMask;
}

/*=============================================================================
函数名    ：CAntImage::setParameters()
功能      ：设置形态学操作的参数
算法实现  ：无
参数说明  ：targetSizeX  
targetSizeY  
xDilateSize  
yDilateSize  
返回值说明：无
其他说明  ：无
=============================================================================*/
void CAntImage::setParameters(int targetSizeX,int targetSizeY, int xDilateSize, int yDilateSize)
{
    m_minTargetSizeX = targetSizeX;
    m_minTargetSizeY = targetSizeY;
    m_x_dilate_size = xDilateSize;
    m_y_dilate_size = yDilateSize;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WIN32
#include "intrinsic.h"
#define restrict
#endif

/*=============================================================================
函数名    ：dilate_erode_HTrans2()
功能      ：膨胀和腐蚀
算法实现  ：形态学膨胀和腐蚀
参数说明  ：pu8Src
m_width
m_height
pu8Dst
halfsizeX
dilate_mask
erode_mask
返回值说明：无
其他说明  ：无
=============================================================================*/
void dilate_erode_HTrans2(u8 * restrict pu8Src, int m_width, int m_height,  
                          u8 * restrict pu8Dst, int halfsizeX, int dilate_mask, int erode_mask)
{
    int i, j, n, k, id;
    int k_last;
    int k_end;
    int i_meet0;
    int i_meetX;
    int u8Cur, u8Prev;
    s64 s64Search, s64Data1, s64Data2, s64Data3, s64Data4;
    u32 u32Dilate, u32Erode, u32Mask, u32Mask1, u32Mask2, u32Mask3, u32Mask4;
    l32 l32BitCount;
    u8 * restrict pu8DstCol;

    int iChangePos[704];
    int iChangeCnt =0;
    int *piChangePos;

    u32Dilate = dilate_mask | (dilate_mask << 8);
    u32Dilate = u32Dilate | (u32Dilate << 16);

    u32Erode = erode_mask | (erode_mask << 8);
    u32Erode = u32Erode | (u32Erode << 16);

    for (i=0; i<m_height; i++)
    {
        pu8DstCol = pu8Dst + i;
#if 1
        piChangePos = iChangePos;
        iChangeCnt = 0;
        u8Prev = erode_mask;    //第一个跳变位置是需要膨胀的mask
        for(j=0;j<m_width;j++)
        {
            u8Cur = pu8Src[j];
            if(u8Cur != u8Prev)
            {
                //记录跳变位置
                *piChangePos++ = j;

            }
            u8Prev = u8Cur;
        }
        *piChangePos =  m_width;//一行结束位置视为一次跳变
        iChangeCnt = piChangePos - iChangePos;

        k =0;
        for(id=0; id<iChangeCnt; id+=2)
        {
            //一个扫描周期分为两个部分：分别用下划线和双下划线标识
            // x x x x x 0 0 0 0 0 ?
            // --------- =========
            //从非零遇到零的时刻，可以确定之前的非零数据复制多少
            //从零遇到非零的时刻，可以确定可以填充多少个零到目标
            i_meet0 = iChangePos[id];
            i_meetX = iChangePos[id+1];

            //从上次目的位置，到i_meet0-halfsizeX位置，复制原始数据
            k_end = i_meet0 - halfsizeX;
            if(i_meet0 == i_meetX) k_end = i_meet0; //非零数据遇到行尾
            for(; k<k_end; k++) 
            {
                //pu8DstCol[0] = pu8Src[k];
                pu8DstCol[0] = erode_mask;//直接使用前景像素点即可，不用复制
                pu8DstCol += m_height;
            }
            //继续填充零到i_meetX + halfsizeX
            k_end = i_meetX + halfsizeX;
            if(k_end > m_width) k_end = m_width;
            for(; k<k_end; k++) 
            {
                pu8DstCol[0] = dilate_mask;
                pu8DstCol += m_height;
            }
        }
        for(; k<m_width; k++) 
        {
            pu8DstCol[0] = erode_mask;
            pu8DstCol += m_height;
        }

#else
        k_last = 0;
        for(j=0;j<m_width;j++)
        {
            //一个扫描周期分为两个部分：分别用下划线和双下划线标识
            // x x x x x 0 0 0 0 0 ?
            // --------- =========
            //从非零遇到零的时刻，可以确定之前的非零数据复制多少
            //从零遇到非零的时刻，可以确定可以填充多少个零到目标

#if 1
            while(pu8Src[j] != dilate_mask && j<m_width) j++;
            i_meet0 = j;
            while(pu8Src[j] == dilate_mask && j<m_width) j++;
            i_meetX = j;
#else
            //使用compeq形成的位图，配合lmbd位检测指令加速寻找目标
            s64Search = _itoll(u32Dilate, u32Dilate);
            l32BitCount = 32;
            while(l32BitCount == 32 && j<m_width)
            {
                s64Data1 = _mem8_const(pu8Src + j);
                s64Data2 = _mem8_const(pu8Src + j + 8);
                s64Data3 = _mem8_const(pu8Src + j + 16);
                s64Data4 = _mem8_const(pu8Src + j + 24);

                u32Mask1 = _dcmpeq4(s64Data1, s64Search);
                u32Mask2 = _dcmpeq4(s64Data2, s64Search);
                u32Mask3 = _dcmpeq4(s64Data3, s64Search);
                u32Mask4 = _dcmpeq4(s64Data4, s64Search);

                u32Mask1 = _pack2(u32Mask2, u32Mask1);
                u32Mask3 = _pack2(u32Mask4, u32Mask3);
                u32Mask3 = _packl4(u32Mask3, u32Mask1);
                u32Mask3 = _bitr(u32Mask3);
                l32BitCount = _lmbd(1, u32Mask3);
                j+= 32;
            }
            j = j - 32 + (l32BitCount);
            if(j>m_width) j=m_width;
            i_meet0 = j;

            s64Search = _itoll(u32Erode, u32Erode);
            l32BitCount = 32;
            while(l32BitCount == 32 && j<m_width)
            {
                s64Data1 = _mem8_const(pu8Src + j);
                s64Data2 = _mem8_const(pu8Src + j + 8);
                s64Data3 = _mem8_const(pu8Src + j + 16);
                s64Data4 = _mem8_const(pu8Src + j + 24);

                u32Mask1 = _dcmpeq4(s64Data1, s64Search);
                u32Mask2 = _dcmpeq4(s64Data2, s64Search);
                u32Mask3 = _dcmpeq4(s64Data3, s64Search);
                u32Mask4 = _dcmpeq4(s64Data4, s64Search);

                u32Mask1 = _pack2(u32Mask2, u32Mask1);
                u32Mask3 = _pack2(u32Mask4, u32Mask3);
                u32Mask3 = _packl4(u32Mask3, u32Mask1);
                u32Mask3 = _bitr(u32Mask3);
                l32BitCount = _lmbd(1, u32Mask3);
                j+= 32;
            }
            j = j - 32 + (l32BitCount);
            if(j>m_width) j=m_width;
            i_meetX = j;
#endif
            //从上次目的位置，到i_meet0-halfsizeX位置，复制原始数据
            k_end = i_meet0 - halfsizeX;
            if(i_meet0 == i_meetX) k_end = i_meet0; //非零数据遇到行尾
            for(k=k_last; k<k_end; k++) 
            {
                //pu8DstCol[0] = pu8Src[k];
                pu8DstCol[0] = erode_mask;//直接使用前景像素点即可，不用复制
                pu8DstCol += m_height;
            }
            //继续填充零到i_meetX + halfsizeX
            k_end = i_meetX + halfsizeX;
            if(k_end > m_width) k_end = m_width;
            for(; k<k_end; k++) 
            {
                pu8DstCol[0] = dilate_mask;
                pu8DstCol += m_height;
            }
            k_last = k;
        }
#endif
        /*
        //如果不是为了保持跟源算法一致，下面的循环可以去掉
        pu8DstCol = pu8Dst + i + (m_width-halfsizeX)*m_height;
        for(j=m_width-halfsizeX;j<m_width; j++)
        {
        pu8DstCol[0] = 0;
        pu8DstCol += m_height;
        }
        */
        pu8Src += m_width;
    }
}

/*=============================================================================
函数名    ：MophOpen_V1()
功能      ：形态学开运算操作
算法实现  ：形态学开运算
参数说明  ：m_pFgMask
m_width
m_height
m_minTargetSizeX  
m_minTargetSizeY  
m_x_dilate_size   
m_y_dilate_size   
m_pTempMask       
返回值说明：无
其他说明  ：无
=============================================================================*/
void MophOpen_V1(u8 *m_pFgMask, int m_width, int m_height, 
                 int m_minTargetSizeX, int m_minTargetSizeY,
                 int m_x_dilate_size, int m_y_dilate_size,
                 u8 *m_pTempMask)
{
    int halfsizeX;
    int halfsizeY;
    if(m_minTargetSizeX > 0 && m_minTargetSizeY > 0) 
    {
        halfsizeX = m_minTargetSizeX/2;
        halfsizeY = m_minTargetSizeY/2;
        dilate_erode_HTrans2(m_pFgMask, m_width, m_height,  m_pTempMask, halfsizeX, 0, 0xFF);
        dilate_erode_HTrans2(m_pTempMask, m_height, m_width,m_pFgMask, halfsizeY, 0, 0xFF);
    }
    
    if(m_x_dilate_size > 0 && m_y_dilate_size > 0)
    {
        halfsizeX = m_x_dilate_size/2;
        halfsizeY = m_y_dilate_size/2;
        dilate_erode_HTrans2(m_pFgMask, m_width, m_height,  m_pTempMask, halfsizeX, 0xFF, 0);
        dilate_erode_HTrans2(m_pTempMask, m_height, m_width,m_pFgMask, halfsizeY,  0xFF, 0);
    }
}


/*=============================================================================
函数名    ：NSET_link()
功能      ：新增一个链接
算法实现  ：无
参数说明  ：pNSET
id1
id2
返回值说明：无
其他说明  ：NSET每个元素有两个特性，id(数组索引)和next_id(数组值)
=============================================================================*/

//新增一个链接
void NSET_link(u32 * pNSET, u32 id1, u32 id2)
{
    int prev_id, next_id;
    int nz_id = 0;
    int other_id;

    if(pNSET[id1]) 
    {
        nz_id = id1;
        other_id = id2;
    }
    else 
        if(pNSET[id2]) 
        {
            nz_id = id2;
            other_id = id1;
        }

        if(nz_id == 0)
        {
            //两个id都是空的，新建一个nset
            pNSET[id1] = id2;
            pNSET[id2] = id1;
        }
        else
        {
            //从非空的nset中遍历寻找是否已经存在跟other_id的链接
            next_id = pNSET[nz_id];
            while(next_id != nz_id)
            {
                if(next_id == other_id) return;

                prev_id = next_id;
                next_id = pNSET[prev_id];
            }

            if(pNSET[other_id])
            {
                //other_id代表一个新的nset,链接两个nset
                //进入otherid代表的nset环路
                pNSET[prev_id] = other_id;
                do
                {
                    prev_id = pNSET[prev_id];
                }while(pNSET[prev_id] != other_id);
                //修改环路输出到nz_id环
                pNSET[prev_id] = nz_id;
            }
            else
            {
                //other_id代表一个单独元素,
                next_id = pNSET[nz_id];
                pNSET[nz_id] = other_id;
                pNSET[other_id] = next_id;
            }
        }
}


/*
对二值化后的图像进行联通区域分析
功能类似matlab的bwlabel
区域分别使用1~N标识在结果数组中

pu8SrcBW (l32Width x l32Height) 二值化后的源图像
u8Fg                            二值源图像的前景色
pu8DstLB                        联通区域标定之后的结果id图，大小也是(l32Width x l32Height)
跟源图像逐点对应，id为0代表背景，为1~N代表联通的前景
返回值N                         联通前景个数。
*/

u32 au32PrevConnID[IMAGE_WIDTH];
u32 au32NSET[IMAGE_WIDTH*IMAGE_HEIGHT/8];

TRunLength atLineRL[IMAGE_WIDTH*IMAGE_HEIGHT/4]; //保存整张图像的RunLength表(由于之前形态学的缘故，此表不可能太大)

/*=============================================================================
函数名    ：findConnectedComponents_V1()
功能      ：对二值化后的图像进行连通区域分析，功能类似matlab的bwlabel，区域分别使用1~N标识在结果数组中
算法实现  ：无
参数说明  ：pu8SrcBW    二值化后的源图像
pu32DstLB   连通区域标定之后的结果id图，大小也是(l32Width x l32Height)
l32Width    图像宽
l32Height   图像高
l32Stride   图像行数据步长
u8Fg        二值源图像的前景色
m_map      
返回值说明：无
其他说明  ：
=============================================================================*/
int findConnectedComponents_V1(u8 *pu8SrcBW, u32 *pu32DstLB, 
                               l32 l32Width, l32 l32Height, l32 l32Stride, u8 u8Fg, int *m_map)
{
    u32 *pu32PrevNB = au32PrevConnID + 1;
    u32 * restrict pu32LBCur;
    int *pl32IDMap = m_map;

    l32 x, y, id, area_id, next_id, overlap, x_st;
    u8 * restrict pu8Src;
    l32 seg_st, seg_ed, seg_id, conn_id;
    l32 l32PrevNBCount, l32GroupID;

    TRunLength * ptRLCur;
    TRunLength * ptRLPrev;
    TRunLength * ptRLNextPrev;
    TRunLength * ptRLTemp;

    s64 s64Data1, s64Data2, s64Data3, s64Data4, s64Search;
    u32 u32Mask1, u32Mask2, u32Mask3, u32Mask4;
    l32 l32BitCount;

    memset(au32NSET, 0, sizeof(au32NSET));

    //自顶向下扫描，综合一行内联通部分，以及垂直方向向下分叉的联通部分
    //为此初步联通区域使用area_id标定，而由于垂直向上分叉引起的已经标定
    //的area的联通，则记录到au8Link数组中，供后级综合之用.
    ptRLCur = atLineRL;
    pu8Src = pu8SrcBW;
    area_id = 0;
    ptRLNextPrev = NULL;
    for(y=0; y<l32Height; y++)
    {
        //记录下一此循环的Prev RunLength信息
        ptRLPrev = ptRLNextPrev;
        ptRLNextPrev = ptRLCur;

        x = 0;
        while(x < l32Width)
        {
            //标志连续的背景，寻找本段seg的起始位置
#if 0
            x_st = x;
            s64Search = _itoll(0xFFFFFFFF, 0xFFFFFFFF);
            l32BitCount = 32;
            while(l32BitCount == 32 && x<l32Width)
            {
                s64Data1 = _mem8_const(pu8Src + x);
                s64Data2 = _mem8_const(pu8Src + x + 8);
                s64Data3 = _mem8_const(pu8Src + x + 16);
                s64Data4 = _mem8_const(pu8Src + x + 24);

                u32Mask1 = _dcmpeq4(s64Data1, s64Search);
                u32Mask2 = _dcmpeq4(s64Data2, s64Search);
                u32Mask3 = _dcmpeq4(s64Data3, s64Search);
                u32Mask4 = _dcmpeq4(s64Data4, s64Search);

                u32Mask1 = _pack2(u32Mask2, u32Mask1);
                u32Mask3 = _pack2(u32Mask4, u32Mask3);
                u32Mask3 = _packl4(u32Mask3, u32Mask1);
                u32Mask3 = _bitr(u32Mask3);
                l32BitCount = _lmbd(1, u32Mask3);
                x+= 32;
            }
            x = x - 32 + (l32BitCount);
            if(x>l32Width) x = l32Width;
            seg_st = x; 
            if(x>=l32Width) break;

            s64Search = _itoll(0, 0);
            l32BitCount = 32;
            while(l32BitCount == 32 && x<l32Width)
            {
                s64Data1 = _mem8_const(pu8Src + x);
                s64Data2 = _mem8_const(pu8Src + x + 8);
                s64Data3 = _mem8_const(pu8Src + x + 16);
                s64Data4 = _mem8_const(pu8Src + x + 24);

                u32Mask1 = _dcmpeq4(s64Data1, s64Search);
                u32Mask2 = _dcmpeq4(s64Data2, s64Search);
                u32Mask3 = _dcmpeq4(s64Data3, s64Search);
                u32Mask4 = _dcmpeq4(s64Data4, s64Search);

                u32Mask1 = _pack2(u32Mask2, u32Mask1);
                u32Mask3 = _pack2(u32Mask4, u32Mask3);
                u32Mask3 = _packl4(u32Mask3, u32Mask1);
                u32Mask3 = _bitr(u32Mask3);
                l32BitCount = _lmbd(1, u32Mask3);
                x+= 32;
            }
            x = x - 32 + (l32BitCount);
            if(x>l32Width) x = l32Width;
            seg_ed = x;

#else
            while((pu8Src[x] == 0) & (x < l32Width))
            {
                x ++;
            }
            if(x == l32Width) break;
            seg_st = x;

            //寻找本端seg的结束位置
            while((pu8Src[x] == u8Fg) & (x < l32Width))
            {
                x++;
            }
            seg_ed = x;
#endif

            //跟前一行seg找邻接关系
            l32PrevNBCount = 0;
            if(y > 0)
            {
                pu32PrevNB[-1] = 0xFFFFFFFF;//此标记用来节省后面while循环中第一个元素的判断
                //检查前一行的RunLength信息判断邻接关系
                ptRLTemp = ptRLPrev;
                while(ptRLTemp->seg_id != 0)
                {
                    overlap = MIN(seg_ed, ptRLTemp->seg_ed) - MAX(seg_st, ptRLTemp->seg_st);
                    if(overlap >= 0) 
                    {
                        if(pu32PrevNB[l32PrevNBCount - 1] != ptRLTemp->seg_id)
                        {
                            pu32PrevNB[l32PrevNBCount] = ptRLTemp->seg_id;
                            l32PrevNBCount ++;
                        }
                    }
                    ptRLTemp ++;
                }
            }

            if(l32PrevNBCount == 0)
            {
                //未在上方找到相邻区域，则此seg对应一个新区域
                area_id ++;
                seg_id = area_id;
            }
            else
            {
                //上方找到连通区域
                seg_id = pu32PrevNB[0];
                for(id = 1; id<l32PrevNBCount; id++)
                {
                    //seg_id跟conn_id存在连接关系，记录下来
                    if(seg_id < sizeof(au32NSET)/sizeof(au32NSET[0]))
                    {
                        NSET_link(au32NSET, seg_id, pu32PrevNB[id]);
                    }
                }
            }

            ptRLCur->seg_id = seg_id;
            ptRLCur->seg_st = seg_st;
            ptRLCur->seg_ed = seg_ed;
            ptRLCur ++;
        }

        //RunLength行结束标志
        ptRLCur->seg_id = 0;
        ptRLCur++;

        pu8Src += l32Stride;
    }

    //根据连接关系，归纳总共组个数(此连接关系数量很少，运算较快)
    //得到一个映射表，输入索引为之前得到的label位图中的area_id,输出为所映射到的组id
    //二者之间是一个多对1的关系。
    //bmp8_write("label_step1.bmp", pDst, 10);

    //NSET数据结构已经把属于相同物体的区域id链接成为一个循环链表
    //只要分别遍历循环链表即可得到映射表
    l32GroupID = 1;
    for(id = 1; id <= area_id; id++)
    {
        pl32IDMap[id] = 0;
    }
    for(id = 1; id <= area_id; id++)
    {
        if(pl32IDMap[id] == 0)
        {
            if(au32NSET[id] == 0)
            {
                //单独一个区域
                pl32IDMap[id] = l32GroupID;
            }
            else
            {
                //多联通区域
                next_id = id;
                do
                {
                    pl32IDMap[next_id] = l32GroupID;
                    next_id = au32NSET[next_id];
                }while(next_id != id);
            }
            l32GroupID ++;
        }
    }
    pl32IDMap[0] = 0;//背景色仍然映射到背景色,我们用0代表背景


#if 1
    //替换RunLength的id为标定过的Label id
    ptRLCur = atLineRL;
    for(y=0; y<l32Height; y++)
    {
        x = 0;
        while(ptRLCur->seg_id > 0)
        {
            id = pl32IDMap[ptRLCur->seg_id];
            ptRLCur->seg_id = id;

            ptRLCur++;
        }
        ptRLCur++;
    }
#else
    //替换id位图为联通元素的组位图,由于替换过程是逐个像素进行的，因此
    //重新特换得到的ID为连续的ID
    pu32LBCur = pu32DstLB;
    ptRLCur = atLineRL;
    for(y=0; y<l32Height; y++)
    {
        x = 0;
        while(ptRLCur->seg_id > 0)
        {
            id = pl32IDMap[ptRLCur->seg_id];

            for(;x<ptRLCur->seg_st;x++) pu32LBCur[x] = 0;
            for(;x<ptRLCur->seg_ed;x++) pu32LBCur[x] = id;

            ptRLCur++;
        }
        for(; x<l32Width; x++) pu32LBCur[x] = 0;

        ptRLCur++;
        pu32LBCur += l32Stride;
    }
#endif
    return l32GroupID;
}


/*=============================================================================
函数名    ：CAntImage::Handle()
功能      ：对前景图像进行形态学操作和连通域分析
算法实现  ：无
参数说明  ：blobs
nBlobs
返回值说明：无
其他说明  ：无
=============================================================================*/
void CAntImage::ExtractBlobByMophCC(u8 * pu8FgMask)                 //输入前景mask图像,经过本函数之后会被膨胀腐蚀修改
{
    int nCCNumber;
    //内部函数依赖此二成员作为输入参数
    m_pFgMask = pu8FgMask;

    //按照设定的参数膨胀腐蚀
    if(0)
    {
        //高架设使用固定经验参数，目标小，因此不腐蚀
        MophOpen_V1(m_pFgMask, m_width, m_height, 
                    0, 0,
                    m_x_dilate_size, m_y_dilate_size,
                    m_pTempMask);
        //
        //目前m_pTempf1没有输出，findBlob2直接使用RunLength信息工作
        nCCNumber = findConnectedComponents_V1(m_pFgMask, (u32*)m_pTempf1,  m_width, m_height, m_width, 255, m_map);

        //findBlobs2函数会根据m_minTargetSizeX, m_minTargetSizeY设定删除过小的目标
        findBlobs2(atLineRL, nCCNumber);
    }
    else
    {
        MophOpen_V1(m_pFgMask, m_width, m_height, 
                    m_minTargetSizeX, m_minTargetSizeY,
                    m_x_dilate_size, m_y_dilate_size,
                    m_pTempMask);
        //目前m_pTempf1没有输出，findBlob2直接使用RunLength信息工作
        nCCNumber = findConnectedComponents_V1(m_pFgMask, (u32*)m_pTempf1,  m_width, m_height, m_width, 255, m_map);
        findBlobs2(atLineRL, nCCNumber);
    }
    return;
}

/*=============================================================================
函数名    ：CAntImage::erodeMask()
功能      ：
算法实现  ：无
参数说明  ：无
返回值说明：无
其他说明  ：
=============================================================================*/
BOOL CAntImage::erodeMask()
{
    int i, j, halfsizeX,halfsizeY, m,n;
    int IsSetToZero;

    halfsizeX = m_minTargetSizeX/2;
    halfsizeY = m_minTargetSizeY/2;

    // first copy foreground mask to temp mask
    memcpy(m_pTempMask, m_pFgMask, sizeof(unsigned char) * m_width * m_height);

    for (i=0; i<m_height; i++)
    {
        for(j=0; j<m_width; j++)
        {
            if(i>=halfsizeY && i<m_height-halfsizeY && j>=halfsizeX && j < m_width-halfsizeX)
            {
                IsSetToZero = 0;
                for(m=-halfsizeY; m<=halfsizeY && !IsSetToZero; m++)
                    for (n=-halfsizeX; n<=halfsizeX &&!IsSetToZero ;n++)
                    {
                        if(m_pFgMask[(i+m)*m_width + (j+n)] == 0)
                        {
                            m_pTempMask[i*m_width+j] = 0; //black
                            IsSetToZero = 1;
                        }
                    } 	
            }
            else
            {
                m_pTempMask[i*m_width+j] = 0; //black
            }
        }
    }

    return TRUE;
}

/*=============================================================================
函数名    ：CAntImage::dilateMask()
功能      ：
算法实现  ：无
参数说明  ：无
返回值说明：
其他说明  ：
=============================================================================*/
BOOL CAntImage::dilateMask()
{
    int i, j, x_halfsize, y_halfsize, m, n;

    // first copy in_mask to out_mask
    memcpy(m_pFgMask,m_pTempMask, sizeof(unsigned char) * m_width * m_height);

    x_halfsize = m_x_dilate_size/2; //horizontal half size
    y_halfsize = m_y_dilate_size/2; //vertical half size 
    for (i=y_halfsize; i<m_height-y_halfsize; i++)
    {
        for(j=x_halfsize;j<m_width-x_halfsize;j++)
        {
            if(m_pTempMask[i*m_width+j])
            {
                for(m=-y_halfsize; m<=y_halfsize; m++)
                    for (n=-x_halfsize; n<=x_halfsize;n++)
                    {
                        m_pFgMask[(i+m)*m_width + (j+n)] = 255; //white
                    }
            }
        }
    }

    return TRUE;
}

/*=============================================================================
函数名    ：CAntImage::findConnectedComponents()
功能      ：
算法实现  ：无
参数说明  ：无
返回值说明：
其他说明  ：
=============================================================================*/
BOOL CAntImage::findConnectedComponents()
{
    unsigned char *mask = (unsigned char*)m_pFgMask;
    int *id = (int*)m_pTempf1;
    int x, y, i, j, k, nextId;
    int ids[4];
    int idOfs[4] = { -m_width-1, -m_width, -m_width+1, -1 };   

    // to start, every id value maps to itself
    nextId = 1;
    for(i=0; i<MAX_CC_IDS; i++) m_map[i] = i;

    // scan first pixel as a special case
    if (*mask) 
        *id = nextId++;
    else 
        *id = 0;
    mask++;
    id++;

    // scan rest of first row as a special case   
    for(x=1; x<m_width; x++)
    {
        if (*mask)
        {
            j = *(id - 1);
            if (j > 0) 
                *id = j;
            else 
                *id = nextId++;
        }
        else *id = 0;

        mask++;
        id++;
    }

    // scan rest of rows
    for(y=1; y<m_height; y++)
    {
        // check first pixel of row as a special case
        if (*mask)
        {
            i = *(id - m_width);
            j = *(id - m_width + 1);

            if (j>i) i = j;
            if (i>0) *id = i;
            else *id = nextId++;
        }
        else *id = 0;
        mask++;
        id++;

        // now check the 'middle' of the row
        for(x=1; x<m_width-1; x++)
        {
            if (*mask)
            {
                j = 0;
                // find the max neighbor
                for(i=0; i<4; i++)
                {
                    k = *(id + idOfs[i]);
                    if(k>=MAX_CC_IDS)
                        ids[i] = 0;
                    else
                        ids[i] = m_map[k];
                    if (ids[i] > j) j = ids[i];
                }

                if (j > 0)
                {
                    for(i=0; i<4; i++)
                    {
                        if (ids[i]==0 || ids[i]==j) continue;
                        for(k=1; k<nextId; k++)
                        {
                            if (m_map[k]==ids[i]) m_map[k] = j;
                        }                  
                    }
                    *id = j;
                }
                else{
                    *id = nextId++;
                }
            }
            else *id = 0;

            mask++;
            id++;
        }

        // finally, we can check the last pixel of the row as a special case
        if (*mask)
        {
            i = *(id - m_width - 1);
            j = *(id - m_width);         
            if (j>i) i = j;

            j = *(id - 1);
            if (j>i) i = j;

            if (i>0) *id = i;
            else *id = nextId++;
        }
        else *id = 0;
        mask++;
        id++;

        if (nextId >= MAX_CC_IDS)
        {
            //printf("Error - not enough connected component ids (%d)\n", MAX_CC_IDS);
            return FALSE;
        }
    }

    // pass 2 - update ids in label image according to equiv map
    id = (int*)m_pTempf1;

    for(i=0; i<m_width*m_height; i++)
    {
        if (*id > 0) *id = m_map[*id];
        id++;
    }

    return TRUE;
}

/*=============================================================================
函数名    ：CAntImage::findBlobs()
功能      ：
算法实现  ：无
参数说明  ：无
返回值说明：
其他说明  ：无
=============================================================================*/
BOOL CAntImage::findBlobs()
{
    int minSize = m_minTargetSizeX * m_minTargetSizeY;
    int *id = (int*)m_pTempf1;
    int x, y, i, j, n;
    AntBlob *tblob;

    m_nBlobs = 0;

    for(i=0; i<MAX_CC_IDS; i++) 
        m_map[i] = -1;

    n = 0;
    for(y=0; y<m_height; y++)
    {
        // blob[X]->user will serve a dual purpose:
        //  1) flag to say whether the blob 'grew' due to the last row
        //  2) index into map array that references this blob
        for(i=0; i<n; i++) 
        {
            m_blobs[i]->user = 0;
        }

        for(x=0; x<m_width; x++)
        {
            j = *id;         
            // is this a blob pixel?
            if (j > 0)
            {
                if (m_map[j] < 0)
                {
                    if (n < MAX_NUM_BLOBS)
                    {
                        // this is a new blob
                        m_map[j] = n;
                        m_blobs[n]->mass = 1;
                        //m_blobs[n]->centroid.x = x;
                        //m_blobs[n]->centroid.y = y;
                        m_blobs[n]->x0 = x;
                        m_blobs[n]->y0 = y;
                        m_blobs[n]->x1 = x;
                        m_blobs[n]->y1 = y;
                        m_blobs[n]->id = j;
                        m_blobs[n]->user = j;
                        n++;
                    }
                }
                else{
                    // this blob already exists
                    m_blobs[m_map[j]]->user = j;
                    j = m_map[j];
                    m_blobs[j]->mass++;
                    //m_blobs[j]->centroid.x += x;
                    //m_blobs[j]->centroid.y += y;
                    if (x > m_blobs[j]->x1) m_blobs[j]->x1 = x;
                    else if (x < m_blobs[j]->x0) m_blobs[j]->x0 = x;

                    if (y > m_blobs[j]->y1) m_blobs[j]->y1 = y;
                    else if (y < m_blobs[j]->y0) m_blobs[j]->y0 = y;               
                }
            }

            id++;
        }      

        // check for finished blobs
        for(i=0; i<n; i++)
        {
            if ((m_blobs[i]->user == 0) && (m_blobs[i]->mass < minSize))
            {
                // kill this blob by moving it to the end of the list and
                //  decrement the list size
                n--;

                // we have to do a proper swap
                tblob = m_blobs[i];
                m_blobs[i] = m_blobs[n];
                m_blobs[n] = tblob;            

                // we also have to update the map            
                m_map[m_blobs[i]->user] = i;

                // since we changed the blob at position <i> we want to reprocess
                //  it next time through the loop
                i--;
            }
        }
    }

    // do some per blob post-processing
    for(i=0; i<n; i++)
    {
        // we accumulated the sum of all (x, y) values
        //  now we need to divide by the mass to get the centroid
        //m_blobs[i]->centroid.x /= m_blobs[i]->mass;
        //m_blobs[i]->centroid.y /= m_blobs[i]->mass;
    }

    m_nBlobs = n;

    return TRUE;
}

/*=============================================================================
函数名    ：CAntImage::findBlobs2()
功能      ：
算法实现  ：无
参数说明  ：ptLineRL
nCCNumber
返回值说明：无
其他说明  ：利用RunLength图生成Blob,逐行完成，同时进行筛选过小Blob，因此不需要使用label图像了
=============================================================================*/
BOOL CAntImage::findBlobs2(TRunLength *ptLineRL, int nCCNumber)
{
    int minSize = m_minTargetSizeX * m_minTargetSizeY;
    int x, y, i, j, n, id;
    TRunLength * ptRLCur;
    int seg_len;
    AntBlob *tblob;

    m_nBlobs = 0;

    for(i = 0; i < nCCNumber; i++) m_map[i] = -1;

    //m_map[i]=j,代表第i个label算法得到的blob被映射到第j个m_blobs结构体
    n = 0;
    ptRLCur = ptLineRL;
    for(y=0; y<m_height; y++)
    {
        for(i=0; i<n; i++) m_blobs[i]->user = 0;

        //第y行数据分析，以连续的前景段为单位
        while(ptRLCur->seg_id > 0)
        {
            id = ptRLCur->seg_id;
            //(ptRLCur->seg_st, ptRLCur->seg_ed)之间就是前景，且其id为id
            //找到其映射到的blob

            if(m_map[id] == -1 && n<MAX_NUM_BLOBS)  //越界判定
            {
                //新建blob
                j = m_map[id] = n;
                n++;

                //初始化此blob
                //m_blobs[j]->centroid.x = 0;
                //m_blobs[j]->centroid.y = 0;
                m_blobs[j]->mass = 0;
                m_blobs[j]->x0 = m_width;
                m_blobs[j]->y0 = m_height;
                m_blobs[j]->x1 = 0;
                m_blobs[j]->y1 = 0;

            }

            if(m_map[id] != -1)      //越界判定
            {
                j = m_map[id];
                m_blobs[j]->user = id;
                //累计R,G,B(灰度在最后折算)
                seg_len = (ptRLCur->seg_ed - ptRLCur->seg_st);
                i = ptRLCur->seg_st * 3;
                /*
                for(x = ptRLCur->seg_st; x<ptRLCur->seg_ed; x++, i+=3)
                {
                    m_blobs[j]->centroid.x += x;
                }
                m_blobs[j]->centroid.y += y*seg_len;
                */
                //更新质量
                m_blobs[j]->mass += seg_len;
                //更新区域
                if ((ptRLCur->seg_ed-1) > m_blobs[j]->x1) m_blobs[j]->x1 = (ptRLCur->seg_ed-1);
                if ((ptRLCur->seg_st)   < m_blobs[j]->x0) m_blobs[j]->x0 = (ptRLCur->seg_st);
                if (y > m_blobs[j]->y1) m_blobs[j]->y1 = y;
                if (y < m_blobs[j]->y0) m_blobs[j]->y0 = y;
            }
            ptRLCur++;
        }

        //本行结束，查找结束的blob，及时删除面积过小者
        for(i=0; i<n; i++)
        {
            if ((m_blobs[i]->user == 0) && (m_blobs[i]->mass < minSize))
            {
                // kill this blob by moving it to the end of the list and decrement the list size
                n--;

                tblob = m_blobs[i];
                m_blobs[i] = m_blobs[n];
                m_blobs[n] = tblob;

                // we also have to update the map            
                m_map[m_blobs[i]->user] = i;

                // since we changed the blob at position <i> we want to reprocess
                //  it next time through the loop
                i--;
            }
        }

        ptRLCur++;
    }

    // do some per blob post-processing
    for(i=0; i<n; i++)
    {
        // we accumulated the sum of all (x, y) values
        //  now we need to divide by the mass to get the centroid
        //m_blobs[i]->centroid.x /= m_blobs[i]->mass;
        //m_blobs[i]->centroid.y /= m_blobs[i]->mass;
    }

    m_nBlobs = n;

    return TRUE;
}

#include <math.h>
