#if 1
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "common.h"
#include "VPMSdk.h"
#include "matrix.h"
#include "wac_ptz.h"
#include "alg_private_bdalg.h"
#include "antWACam.h"
#include "get_roi.h"
#include "show_img.h"

#ifdef WIN32
#include "intrinsic.h"
#include "opencv2/legacy/legacy.hpp"
#include "opencv2/opencv.hpp"
using namespace cv;
#endif

extern GlobalSettings g_GblobalSettings;

static u32 g_au32IIimage[(IMAGE_WIDTH + 1)*(IMAGE_HEIGHT + 1)];


extern "C" u64 UmpGetHTime(void);

l32 KLTTrackingPoint(u8 *pu8Pre, u8 *pu8Cur, l32 l32Width, l32 l32Height,
                     f32 f32PrevX, f32 f32PrevY, f32 &f32CurX, f32 &f32CurY,
                     l32 l32WndW, l32 l32WndH, f32 *pf32Error);

l32 KLTTrackingPointWithWeight(u8 *pu8Pre, u8 *pu8Cur, u8 *pu8PreWeight, l32 l32Width, l32 l32Height,
                               f32 f32PrevX, f32 f32PrevY, f32 &f32CurX, f32 &f32CurY,
                               l32 l32WndW, l32 l32WndH, f32 *pf32Error);
l32 KLTTrackingPointWithWeight2(u8 *pu8Pre, u8 *pu8Cur, u8 *pu8PreWeight, l32 l32Width, l32 l32Height,
								f32 f32PrevX, f32 f32PrevY, f32 &f32CurX, f32 &f32CurY,
								l32 l32WndW, l32 l32WndH, f32 *pf32Error);

void RectShrink(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs);
static void RectExpand(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs);

void sobel_sum_hv(u8 *pu8Src, l32 l32Width, l32 l32Height, l32 l32Stride,
                  l32 &l32SobelPowerV, l32 &l32SobelPowerH);
int isInROI(int x, int y,int *iRoiNum);
static BOOL IsExistFg(TRect *ptObjRc, u8 *pu8FgMask);

// 高度扩展ok
static void RectExpandH(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs, l32 l32MaxHeight)
{
    int ys, ye;
    int consCnt;

    memcpy((void *)pNewRect, (void *)pOldRect, sizeof(TRect));

    // 	l32 al32GrayMax[IMAGE_HEIGHT];
    // 	ys = pOldRect->l32Top - l32MaxHeight;
    // 	ys = ys < 0 ? 0 : ys;
    // 	ye = pOldRect->l32Top + pOldRect->l32Height + l32MaxHeight;
    // 	ye = ye > (IMAGE_HEIGHT-1)?(IMAGE_HEIGHT-1):ye;
    // 	for(int y = ys; y < ye; y++)
    // 	{
    // 		u8 u8GrayMax = 0;
    // 		for(int x = pOldRect->l32Left; x < pOldRect->l32Left + pOldRect->l32Width; x++)
    // 		{
    // 			if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
    // 			{
    // 				al32GrayMax[y - ys] = pu8DiffAbs[y * IMAGE_WIDTH + x];
    // 			}
    // 		}
    // 	}

    // 计算向上扩出的h/2每一行的最值投影
    ys = pOldRect->l32Top;
    ye = pOldRect->l32Top - l32MaxHeight;
    ye = ye < 0 ? 0 : ye;
    consCnt = 0;

    for(int y=ys; y>=ye; y--)
    {
        u8 u8GrayMax = 0;
        for(int x = pOldRect->l32Left; x < pOldRect->l32Left + pOldRect->l32Width; x++)
        {
            if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
            {
                u8GrayMax = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
        if(u8GrayMax > 50)
        {
            pNewRect->l32Top--;
            pNewRect->l32Height++;
            consCnt = 0;
        }
        else
        {
            consCnt++;
        }
        if(consCnt > 5)
        {
            break;
        }
    }
    //  	if(pNewRect->l32Height >= l32MaxHeight)
    //  	{
    //  		return;
    //  	}

    ys = pNewRect->l32Top + pNewRect->l32Height - 1;
    ye = ys + l32MaxHeight;
    ye = ye > (IMAGE_HEIGHT-1)?(IMAGE_HEIGHT-1):ye;
    consCnt = 0;
    for(int y=ys; y<ye; y++)
    {
        u8 u8GrayMax = 0;
        for(int x = pOldRect->l32Left; x < pOldRect->l32Left + pOldRect->l32Width; x++)
        {
            if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
            {
                u8GrayMax = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
        if(u8GrayMax > 50)
        {
            pNewRect->l32Height++;
            consCnt = 0;
        }
        else
        {
            consCnt++;
        }

        if(consCnt > 5)
            break;
    }
}

// 高度缩小ok
static void RectShrinkH(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs, l32 l32MinHeight)
{
    l32 al32Sum[IMAGE_HEIGHT] = {-1};
    int ys, ye;

    memcpy((void *)pNewRect, (void *)pOldRect, sizeof(TRect));
    // 计算每一行的最值投影
    for(int y = pOldRect->l32Top; y < pOldRect->l32Top + pOldRect->l32Height - 1; y++)
    {
        for(int x = pOldRect->l32Left; x < pOldRect->l32Left + pOldRect->l32Width; x++)
        {
            if(al32Sum[y - pOldRect->l32Top] < pu8DiffAbs[y * IMAGE_WIDTH + x])
            {
                al32Sum[y- pOldRect->l32Top] = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
    }

    // 从上往下收缩
    for(int y = 0; y < pOldRect->l32Height - l32MinHeight; y++)
    {
        if(al32Sum[y] < 50)
        {
            pNewRect->l32Top++;
            pNewRect->l32Height--;
        }
        else
        {
            break;
        }
    }
    // 	if(pNewRect->l32Height <= l32MinHeight)
    // 	{
    // 		return;
    // 	}

    // 从下往上收缩
    for(int y = pOldRect->l32Height - 1; y > l32MinHeight - 1; y--)
    {
        if(al32Sum[y] < 50)
        {
            pNewRect->l32Height--;
        }
        else
        {
            break;
        }
    }
}

// 宽度扩展ok
static void RectExpandW(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs, l32 l32MaxWidth)
{
    int xs, xe;
    memcpy((void *)pNewRect, (void *)pOldRect, sizeof(TRect));
    // 	if(pOldRect->l32Width >= l32MaxWidth)
    // 		return;

    // 计算向左扩出的w每一列的最值投影
    xs = pOldRect->l32Left;
    xe = pOldRect->l32Left - l32MaxWidth/*pOldRect->l32Width)*/;
    xe = xe < 0 ? 0 : xe;
    for(int x=xs; x>=xe; x--)
    {
        u8 u8GrayMax = 0;
        for(int y = pOldRect->l32Top; y < pOldRect->l32Top + pOldRect->l32Height; y++)
        {
            if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
            {
                u8GrayMax = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
        if(u8GrayMax > 50)
        {
            pNewRect->l32Left--;
            pNewRect->l32Width++;
        }
        else
        {
            break;
        }
    }
    // 	if(pNewRect->l32Width >= l32MaxWidth)
    // 	{
    // 		return;
    // 	}

    // 按照最大的车的宽度，右边还有多少剩余的宽度余量
    xs = pNewRect->l32Left + pNewRect->l32Width - 1;
    xe = xs + l32MaxWidth - 1;
    xe = xe > (IMAGE_WIDTH - 1) ? (IMAGE_WIDTH - 1) : xe;
    for(int x=xs; x<xe; x++)
    {
        u8 u8GrayMax = 0;
        for(int y = pOldRect->l32Top; y < pOldRect->l32Top + pOldRect->l32Height; y++)
        {
            if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
            {
                u8GrayMax = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
        if(u8GrayMax > 50)
        {
            pNewRect->l32Width++;
        }
        else
        {
            break;
        }
    }
}

// 宽度缩小
static void RectShrinkW(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs, l32 l32MinWidth)
{
    int x, y;
    l32 al32Max[IMAGE_WIDTH] = {-1};
    memcpy((void *)pNewRect, (u8 *)pOldRect, sizeof(TRect));

    if(pOldRect->l32Width <= l32MinWidth)
    {
        return;
    }

    // 每一列的最值投影
    for(x=pOldRect->l32Left; x<pOldRect->l32Left + pOldRect->l32Width-1; x++)
    {
        for(y = pOldRect->l32Top; y < pOldRect->l32Top + pOldRect->l32Height-1; y++)
        {
            if(al32Max[x-pOldRect->l32Left] < pu8DiffAbs[y * IMAGE_WIDTH + x])
            {
                al32Max[x-pOldRect->l32Left] = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
    }

    for(x=0; x<pOldRect[0].l32Width - l32MinWidth; x++)
    {
        if(al32Max[x] < 50)
        {
            pNewRect->l32Left++;
            pNewRect->l32Width--;
        }
        else
        {
            break;
        }
    }
    if(pNewRect->l32Width <= l32MinWidth)
    {
        return;
    }

    for(x=pNewRect->l32Width - 1; x >= l32MinWidth; x--)
    {
        if(al32Max[x] < 50)
        {
            pNewRect->l32Width--;
        }
        else
        {
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//简单kalman跟踪器
KalmanFilter8x3::KalmanFilter8x3()
{
    int i,j;
    stateVec = stateVec_data;
    vecHx = vecHx_data;
    //临时空间
    Mtx_Buff8x8x8 = alloc_DoubleMatrix(ANT_STATE_DIM*8,ANT_STATE_DIM);
    //一次性分配好 (8+3+3)x8的矩阵，通过偏移得到后继矩阵
    stateCovMtx = alloc_DoubleMatrix(ANT_STATE_DIM + ANT_SPACE_DIM + ANT_SPACE_DIM,ANT_STATE_DIM);
    invS = stateCovMtx + ANT_STATE_DIM;
    S = invS + ANT_SPACE_DIM;

    F_Mtx=alloc_DoubleMatrix(ANT_STATE_DIM,ANT_STATE_DIM);
    Q_Mtx=alloc_DoubleMatrix(ANT_STATE_DIM,ANT_STATE_DIM);
    H_Mtx=alloc_DoubleMatrix(ANT_SPACE_DIM,ANT_STATE_DIM);
    R_Mtx=alloc_DoubleMatrix(ANT_SPACE_DIM,ANT_SPACE_DIM);
    P_init_Mtx = alloc_DoubleMatrix(ANT_STATE_DIM,ANT_STATE_DIM);
    /* initialize F matrix to
    //取虚拟T=1
    //x, dx, ax,  y, dy, ay, m, dm
    x 1,  T,
    dx 0,  1,  T
    ax 0,  0,  1,
    y 0,  0,  0,  1,  T
    dy 0,  0,  0,  0,  1,  T
    ay 0,  0,  0,  0,  0,  1,
    m 0,  0,  0,  0,  0,  0, 1, T
    dm 0,  0,  0,  0,  0,  0, 0, 1    //m倾向于不变

    }; //the "F" matrix in target motion model */
    for(i=0;i<ANT_STATE_DIM;i++)
        for(j=0;j<ANT_STATE_DIM;j++)
            if(i==j)
                F_Mtx[i][j]=1.0f;
            else
                F_Mtx[i][j]=0.0;
    F_Mtx[0][1]=1;
    F_Mtx[1][2]=1;
    F_Mtx[3][4]=1;
    F_Mtx[4][5]=1;
    F_Mtx[6][7]=1;
    /* initialize H matrix to be
    H = [1 0 0 0 0 0 0 0 ;
    0 0 0 1 0 0 0 0 ;
    0 0 0 0 0 0 1 0];
    观察矩阵，从状态矩阵中取x,y,m
    */
    for(i=0;i<ANT_SPACE_DIM;i++)
        for(j=0;j<ANT_STATE_DIM;j++)
            H_Mtx[i][j]=0.0;
    H_Mtx[0][0]=1.0f;
    H_Mtx[1][3]=1.0f;
    H_Mtx[2][6]=1.0f;
    /*
    过程噪声协方差矩阵：对角阵
    */
#define TRANS_ERROR (1e-8)
    for(i=0;i<ANT_STATE_DIM;i++)
        for(j=0;j<ANT_STATE_DIM;j++)
            if(i==j)
                Q_Mtx[i][j]=TRANS_ERROR;
            else
                Q_Mtx[i][j]=0.0;
    /* initialize R matrix to be
    R = [sigmaObsSqr	0				0					;
    0				sigmaObsSqr		0					;
    0				0				m_massObsVarianceFactor	];
    assuming
    1) the errors in x and y position are independent
    2) the std of error in x (the same for y axis) is 1 pixel, can be adjusted according to other error models
    观察误差矩阵R
    */
    for(i=0;i<ANT_SPACE_DIM;i++)
        for(j=0;j<ANT_SPACE_DIM;j++)
            R_Mtx[i][j] = 0.0;
    //质心位置的观测误差比较小，mass的误差比较大
    //大尺寸的物体移动能力强，其观测值误差较小
    R_Mtx[0][0] = 1e-8;
    R_Mtx[1][1] = 1e-8;
    R_Mtx[2][2] = 1e1;


    /* The state covariance matix using a two-point differencing track start method
    状态误差协方差矩阵初始值
    */
    for(i=0;i<ANT_STATE_DIM;i++)
        for(j=0;j<ANT_STATE_DIM;j++)
            if (i==j)
                P_init_Mtx[i][j]=1;
            else
                P_init_Mtx[i][j]=0.0;
}

KalmanFilter8x3::~KalmanFilter8x3()
{
    free_DoubleMatrix(Mtx_Buff8x8x8);
    free_DoubleMatrix(stateCovMtx);
}
void KalmanFilter8x3::correct_xy(float x, float y)
{
    stateVec[0] = x;
    stateVec[3] = y;
}
void KalmanFilter8x3::reset(double centroidx, double centroidy, 
                            double centroidxv = 0, double centroidyv = 0,
                            double mass = 0)
{
    //state vector from 2-point differencing
    stateVec[0]=centroidx; //x position
    stateVec[1]=centroidxv; //x velocity
    stateVec[2]=0;
    stateVec[3]=centroidy; //y position
    stateVec[4]=centroidyv; //y velocity
    stateVec[5]=0;
    stateVec[6]=(mass);
    stateVec[7]=0;
    //state covariance matrix is set to P_init_Mtx, using two-point differencing track start method
    CopyDoubleMatrix(ANT_STATE_DIM, ANT_STATE_DIM, P_init_Mtx, &(stateCovMtx));
}

void KalmanFilter8x3::propagateState(int n, double* xkn, double* ykn, double* xk, double* yk)
{
    int i;
    double projState_data[ANT_STATE_DIM];
    double tmpVec_data[ANT_STATE_DIM];
    DoubleVector projState = projState_data;
    DoubleVector tmpVec = tmpVec_data;

    CopyDoubleVector(ANT_STATE_DIM, stateVec, &tmpVec);
    for (i=0;i<n;i++)
    {
        //propagate by 1 frame, save the propagated state in projState
        DoubleMatrixVectorProduct(ANT_STATE_DIM, ANT_STATE_DIM, F_Mtx, tmpVec, &projState);
        //copy projState to tmpVec
        CopyDoubleVector(ANT_STATE_DIM, projState, &tmpVec);
    }
    *xkn=projState[0];
    *ykn=projState[3];
    *xk=stateVec[0];
    *yk=stateVec[3];
}
void KalmanFilter8x3::predictState(float &centroidx, float &centroidy, float &mass)
{
    DoubleMatrix tmpMtxP1; //same dimension as P
    DoubleMatrix tmpMtxP2;
    DoubleMatrix tmpMtxS1; //same dimension as S
    DoubleMatrix tmpMtxH1;
    double stateVec_copy[ANT_STATE_DIM];
    DoubleVector tmpstateVec = (DoubleVector)stateVec_copy;
    tmpMtxP1 = Mtx_Buff8x8x8;
    tmpMtxP2 = tmpMtxP1 + ANT_STATE_DIM;
    tmpMtxS1 = tmpMtxP2 + ANT_STATE_DIM;
    tmpMtxH1 = tmpMtxS1 + ANT_SPACE_DIM;
    //prediction stage
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //propagation of the state vector
    //x(k|k=1)=F*x(k-1|k-1)
    CopyDoubleVector(8, stateVec, &(tmpstateVec));
    DoubleMatrixVectorProduct(ANT_STATE_DIM, ANT_STATE_DIM, F_Mtx, tmpstateVec, &(stateVec));
    //propagation of the covariance matrix
    //Tmp1=F*P(k-1|k-1)
    DoubleMatrixMatrixProduct(ANT_STATE_DIM, ANT_STATE_DIM, ANT_STATE_DIM, F_Mtx, stateCovMtx, &tmpMtxP1);
    //Tmp2=Tmp1*F'=F*P(k-1|k-1)*F'
    DoubleMatrixMatrixTransposeProduct(ANT_STATE_DIM, ANT_STATE_DIM, ANT_STATE_DIM, tmpMtxP1, F_Mtx, &tmpMtxP2);
    //P(k|k-1)=Tmp2+Q=F*P(k-1|k-1)*F'+Q
    DoubleMatrixSum(ANT_STATE_DIM, ANT_STATE_DIM, tmpMtxP2, Q_Mtx, &(stateCovMtx));
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //calculation of the covariance matrix S of the innovation v=z-H*x
    //TmpH1=H*P(k|k-1)
    DoubleMatrixMatrixProduct(ANT_SPACE_DIM, ANT_STATE_DIM, ANT_STATE_DIM, H_Mtx, stateCovMtx, &tmpMtxH1);
    //TmpS1=TmpH1*H'=H*P(k|k-1)*H'
    DoubleMatrixMatrixTransposeProduct(ANT_SPACE_DIM, ANT_STATE_DIM, ANT_SPACE_DIM, tmpMtxH1, H_Mtx, &tmpMtxS1);
    //S(k)=TmpS1+R=H*P(k|k-1)*H'+R
    DoubleMatrixSum(ANT_SPACE_DIM, ANT_SPACE_DIM, tmpMtxS1, R_Mtx, &(S));
    //invS=inverse(S)
    //SymmetricDoubleMatrixInverse(ANT_SPACE_DIM, newState->S, &(newState->invS));
    SymmetricDoubleMatrixInverse3x3(S, &(invS));
    //calculate H*x 给出预测的观察值vecHx[0],vecHx[1],vecHx[2]
    DoubleMatrixVectorProduct(ANT_SPACE_DIM, ANT_STATE_DIM, H_Mtx, stateVec, &(vecHx));

    centroidx = vecHx[0];
    centroidy = vecHx[1];
    mass = vecHx[2];
}
void KalmanFilter8x3::updateState(float &centroidx, float &centroidy, float &mass)
{
    double tmp_z_data[ANT_SPACE_DIM];
    double vecKv_data[ANT_STATE_DIM];

    DoubleMatrix tmpMtxP1 = Mtx_Buff8x8x8; //same dimension as P
    DoubleMatrix tmpMtxP2 = tmpMtxP1 + ANT_STATE_DIM;
    DoubleMatrix tmpMtxP3 = tmpMtxP2 + ANT_STATE_DIM;
    DoubleMatrix K        = tmpMtxP3 + ANT_STATE_DIM;
    DoubleMatrix tmpK     = K + ANT_STATE_DIM;
    DoubleMatrix tmpSpaceMtx = tmpK + ANT_STATE_DIM;
    DoubleVector tmp_z      = tmp_z_data;
    DoubleVector vecKv      = vecKv_data;

    //get current innovation from obs v=z-H*x
    tmp_z[0] = centroidx - vecHx[0];
    tmp_z[1] = centroidy - vecHx[1];
    tmp_z[2] = mass - vecHx[2];

    //TmpK=P(k|k-1)*H'
    DoubleMatrixMatrixTransposeProduct(ANT_STATE_DIM, ANT_STATE_DIM, ANT_SPACE_DIM, stateCovMtx, H_Mtx, &tmpK);
    //the Kalman gain K=tmpK*invS=P(k|k-1)*H'*invS;
    DoubleMatrixMatrixProduct(ANT_STATE_DIM, ANT_SPACE_DIM, ANT_SPACE_DIM, tmpK, invS, &K);
    //calculate K*v
    DoubleMatrixVectorProduct(ANT_STATE_DIM, ANT_SPACE_DIM, K, tmp_z, &vecKv);
    //x(k|k)=x(k|k-1)+K*v, this is the updated state vector
    AddToDoubleVector(ANT_STATE_DIM, vecKv, &(stateVec));
    //tmpP1=K*H
    DoubleMatrixMatrixProduct(ANT_STATE_DIM, ANT_SPACE_DIM, ANT_STATE_DIM, K, H_Mtx, &tmpMtxP1);
    //tmpMtxP2=tmpMtxP1*P(k|k-1)=K*H*P(k|k-1)
    DoubleMatrixMatrixProduct(ANT_STATE_DIM, ANT_STATE_DIM, ANT_STATE_DIM, tmpMtxP1, stateCovMtx, &tmpMtxP2);
    //tmpMtxP1=P(k|k-1)-K*H*P(k|k-1)
    DoubleMatrixDifference(ANT_STATE_DIM, ANT_STATE_DIM, stateCovMtx, tmpMtxP2, &tmpMtxP1);
    //P(k|k)=tmpMtxP1=P(k|k-1)-K*H*P(k|k-1)
    CopyDoubleMatrix(ANT_STATE_DIM, ANT_STATE_DIM, tmpMtxP1, &(stateCovMtx));

    //返回更新后的可观察到状态
    centroidx = stateVec[0];
    centroidy = stateVec[3];
    mass = stateVec[6];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


void findShoulderBorder(u8 *pu8DiffAbs, l32 l32Width, l32 l32Height, 
                        int x0, int y0, int x1, int y1,
                        int &leftShoulderOffset,int &rightShoulderOffset)
{
    //寻找上半身左右边界位置,如果是真实人的话，这个位置不会太远
    int x, y;
    int power, max_power;
    int cw = 2*(x1 - x0);

    //寻找最大能量
    max_power = 0;
    for(x = x0; x < x1; x++)
    {
        u8 *pu8Cur = pu8DiffAbs + y0 * l32Width + x;
        power = 0;
        for(y = y0; y < y1; y++)
        {
            power += pu8Cur[0];
            pu8Cur+=l32Width;
        }
        if(power > max_power) max_power = power;
    }

    int cnt = 0;
    for(x = x0; x > x0 - cw && x>0; x--)
    {
        u8 *pu8Cur = pu8DiffAbs + y0 * l32Width + x;
        power = 0;
        for(y = y0; y < y1; y++)
        {
            power += pu8Cur[0];
            pu8Cur+=l32Width;
        }
        if(power < max_power/2) cnt ++;
        else cnt = 0;

        if(cnt > 1) break;
    }
    leftShoulderOffset = x0-x;

    cnt = 0;
    for(x = x1; x < x1 + cw && x<l32Width; x++)
    {
        u8 *pu8Cur = pu8DiffAbs + y0 * l32Width + x;
        power = 0;
        for(y = y0; y < y1; y++)
        {
            power += pu8Cur[0];
            pu8Cur+=l32Width;
        }
        if(power < max_power/2) cnt ++;
        else cnt = 0;

        if(cnt > 1) break;
    }
    rightShoulderOffset = x-x1;

}
void sobel_sum_HeadShoulder(u8 *pu8Src, l32 l32Width, l32 l32Height, l32 l32Stride,
                            l32 &l32PowerLS, l32 &l32PowerRS, l32 &l32SepPos)
{
    l32 l32OffsetU, l32OffsetD, l32OffsetL, l32OffsetR;
    l32 l32I, l32J;
    l32 l32TempLS, l32TempRS, l32MinPower, l32Max;
    l32 l32PowerSumLS, l32PowerSumRS;

    static l32 al32PowerLS[IMAGE_WIDTH];
    static l32 al32PowerRS[IMAGE_WIDTH];
    l32 l32SobelPowerLS, l32SobelPowerRS;

    for(l32J = 0; l32J < l32Width; l32J++)
    {
        al32PowerLS[l32J] = al32PowerRS[l32J] = 0;
    }
    l32PowerSumLS = l32PowerSumRS = 0;

    u8 as8TempImage[704*396]={255};
    u8 *pu8Cur = pu8Src + l32Stride;
    for(l32I = 1; l32I < l32Height - 1; l32I++, pu8Cur+=l32Stride)
    {
        for(l32J = 1; l32J < l32Width - 1; l32J++)
        {
            l32TempLS = (pu8Cur[-l32Stride + (l32J)] + 2*pu8Cur[-l32Stride + (l32J-1)] + pu8Cur[l32J-1])
                -
                (pu8Cur[l32J+1] + 2*pu8Cur[l32Stride + (l32J+1)] + pu8Cur[l32Stride + (l32J)]);

            l32TempRS = (pu8Cur[-l32Stride + (l32J)] + 2*pu8Cur[-l32Stride + (l32J+1)] + pu8Cur[l32J+1])
                -
                (pu8Cur[l32J-1] + 2*pu8Cur[l32Stride + (l32J-1)] + pu8Cur[l32Stride + (l32J)]);

            if(l32TempLS < 0) l32TempLS = -l32TempLS;
            if(l32TempRS < 0) l32TempRS = -l32TempRS;

            l32TempLS = MIN(l32TempLS, 255);
            l32TempRS = MIN(l32TempRS, 255);

            as8TempImage[l32I * l32Width*3 + l32J] = pu8Cur[l32J];
            as8TempImage[l32I * l32Width*3 + l32J+l32Width] = l32TempLS>255?255:l32TempLS;
            as8TempImage[l32I * l32Width*3 + l32J+l32Width*2] = l32TempRS>255?255:l32TempRS;


            al32PowerLS[l32J] += l32TempLS;
            al32PowerRS[l32J] += l32TempRS;

            l32PowerSumLS += l32TempLS;
            l32PowerSumRS += l32TempRS;
        }
    }
    ShowY8InOpenCV("src", as8TempImage, l32Width*3, l32Height, l32Width*3, 1);



    l32SobelPowerRS = l32SobelPowerLS = 0;
    l32I = l32Max = 0;
    for(l32J = 0; l32J < l32Width; l32J++)
    {
        l32SobelPowerLS += al32PowerLS[l32J];
        l32SobelPowerRS += al32PowerRS[l32J];

        //左右肩分别在各自区域内的区分性,明显性
        l32 l32PowerDiscrimLS = l32SobelPowerLS - (l32SobelPowerRS);
        l32 l32PowerDiscrimRS = (l32PowerSumRS - l32SobelPowerRS) - (l32PowerSumLS - l32SobelPowerLS);


        l32MinPower = MIN(l32PowerDiscrimLS, l32PowerDiscrimRS);
        if(l32Max < l32MinPower)
        {
            l32Max = l32MinPower;
            l32I = l32J;
            l32PowerLS = l32PowerDiscrimLS;
            l32PowerRS = l32PowerDiscrimRS;
        }
    }
    l32SepPos = l32I;


    l32PowerLS = l32PowerSumLS;
    l32PowerRS = l32PowerSumRS;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//光流跟踪器


TrackingTarget::TrackingTarget()
{
    m_l32Life = -1;
    m_bIsMissing = TRUE;
    m_targetid = -1;
    m_cx_init = 0;
    m_cy_init = 0;
	m_last_cw = IMAGE_WIDTH;
	m_last_ch = IMAGE_HEIGHT;
};
TrackingTarget::~TrackingTarget()
{	
};

//使用二次函数，最小二乘法拟和运动模型
void TrackingTarget::fitMovementModel()
{
    m_QLLSx.fitModel();
    m_QLLSy.fitModel();
}
//预测n帧之后的位置和大小---使用最近拟和的模型
void TrackingTarget::PredictLocation(int n, 
                                     double* xkn, double* ykn, 
                                     l32* tw, l32* th)
{
//    propagateState(n, xkn, ykn, xk, yk);
//     *tw = g_GblobalSettings.VanishInfo.GetW(m_Wnorm, *ykn);
//     *th = g_GblobalSettings.VanishInfo.GetH(m_Hnorm, *ykn);
    float x = m_QLLSx.predict_acc_speed(n);
    float y = m_QLLSy.predict_acc_speed(n);

	if(xkn) *xkn = x;
	if(ykn) *ykn = y;

	if(tw || th)
	{
		float w, h;
		h = g_GblobalSettings.VanishInfo.GetH_byCentY(m_Hnorm, y);
		w = g_GblobalSettings.VanishInfo.GetW(m_Wnorm,  y + h/2);
		if(tw) *tw = w;
		if(th) *th = h;
	}
}
void TrackingTarget::FindRectByMeanshift(float &fmx, float &fmy, float fmw, float fmh, 
                                         u8 *pu8Mask, l32 l32Width, l32 l32Height)
{
    int x, y;
    int y0, y1;
    int x0, x1;
    //mean-shift方式可以避免形态学，二值化等操作，直接从可信度图像上进行跟踪
    float f_dist_conv_th = 0.9f;
    float dist_diff;
    float fmx_org = fmx;
    float fmy_org = fmy;
    u8 * pu8Temp;

    // 计算搜索的边界
    int left = fmx - fmw;
    int right = fmx + fmw;
    int top = fmy - fmh;
    int bottom = fmy + fmh;
    do
    {
        x0 = fmx-fmw/2; if(x0<left) x0 = left; //if(x0<0) x0 = 0;
        x1 = fmx+fmw/2; if(x1>right) x1 = right;//if(x1>l32Width) x1=l32Width;

        y0 = fmy-fmh/2; if(y0<top) y0 = top;//if(y0<0) y0 = 0;
        y1 = fmy+fmh/2; if(y1>bottom) y1=bottom;//if(y1>l32Height) y1=l32Height;

        float sum_wx=0, sum_wy=0;
        float sum_w=0;
        //使用flat kernel
        pu8Temp = pu8Mask + y0 * l32Width;
        for(y=y0; y<y1; y++)
        {
            for(x=x0;x<x1;x++)
            {
                if(pu8Temp[x])
                {
                    sum_wx += x * pu8Temp[x];
                    sum_wy += y * pu8Temp[x];
                    sum_w += pu8Temp[x];
                }
            }
            pu8Temp += l32Width;
        }
        if(sum_w < 1)
        {
            //跟踪失败
            fmx = fmx_org;
            fmy = fmy_org;
            break;
        }
        float newx = sum_wx/sum_w;
        float newy = sum_wy/sum_w;
        dist_diff = (newx - fmx)*(newx - fmx) + (newy - fmy)*(newy - fmy);
        fmx = newx;
        fmy = newy;
    }while(dist_diff > f_dist_conv_th*f_dist_conv_th);
}
BOOL TrackingTarget::IsCurTargetValid()
{
    return m_l32Life >= 0 && (!m_bIsMissing);
}

BOOL TrackingTarget::reset(TRect *ptRect, u8 *pu8DiffAbs,u8 *pu8MaskGmm, u8 *pu8Yuv, BOOL bSemiAuto)
{
    float fold_cx = m_cx;
    float fold_cy = m_cy;
    float fold_cw = m_cw;
    float fold_ch = m_ch;
    float fold_DiffAvg = m_f32DiffAvg;
    GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;
    GSTrackSettingInfo      &TrackSettingInfo = g_GblobalSettings.TrackSettingInfo;

    u8 *pu8RoiMask = g_GblobalSettings.ROIInfo.apROIMaskImage[3];
    l32 l32SobelPowerV, l32SobelPowerH;
    u8 *pu8Src;
    TRect atRect[3];

    m_cx = ptRect->l32Left + ptRect->l32Width / 2;
    m_cy = ptRect->l32Top + ptRect->l32Height / 2;
    m_cw = ptRect->l32Width;
    m_ch = ptRect->l32Height;

//计算在输入跟踪框位置时，人的模板大小
    float fhuman_w = VanishInfo.GetW(VanishInfo.f32HumanNormW, m_cy + m_ch/2);//VanishInfo.f32MinW * (m_cy + m_ch/2 - VanishInfo.f32Y0)/(VanishInfo.f32TargetY - VanishInfo.f32Y0);
    float fhuman_h = VanishInfo.GetH(VanishInfo.f32HumanNormH, m_cy + m_ch/2);//VanishInfo.f32MinH * (m_cy + m_ch/2 - VanishInfo.f32Y0)/(VanishInfo.f32TargetY - VanishInfo.f32Y0);
   
//	if(bSemiAuto != 2)

	//仅仅跟踪人时，修改外部传入的尺寸
    if((TrackSettingInfo.m_WorkingMode % 10) == WMODE_PERSON)
    {
        float fhuman_w_lo = fhuman_w;
        float fhuman_h_lo = fhuman_h;
        float fhuman_w_hi = fhuman_w * 1.5;
        float fhuman_h_hi = fhuman_h * 1.5;
        if(m_cw < fhuman_w_lo) m_cw = fhuman_w_lo;
        if(m_ch < fhuman_h_lo) m_ch = fhuman_h_lo;
        if(m_cw > fhuman_w_hi) m_cw = fhuman_w_hi;
        if(m_ch > fhuman_h_hi) m_ch = fhuman_h_hi;

		if(ptRect->l32Width > (l32)(fhuman_w_hi * 1.5))
		{
			goto RESET_FAIL;
		}
    }

    atRect[0].l32Left = m_cx - m_cw / 2;
    atRect[0].l32Top = m_cy - m_ch / 2;
    atRect[0].l32Width = m_cw;
    atRect[0].l32Height = m_ch;

    l32 l32InputCenterX = m_cx;
    l32 l32InputCenterY = m_cy; 
    
    // 在两倍宽高的区域去纠正
	if(bSemiAuto != 2)
	{
		FindRectByMeanshift(m_cx, m_cy, m_cw, m_ch, pu8DiffAbs, IMAGE_WIDTH, IMAGE_HEIGHT);
	}

	// meanshift有时会找到不正确的位置
	if(m_cx < 0 || m_cy < 0 || m_cw < 0 || m_ch < 0)
		goto RESET_FAIL;

    if(!bSemiAuto)
    {
        //仅仅跟踪人时，利用背景减除图像,寻找物体左右真实边缘(也就是左右肩膀,从而估计物体宽度是否合理)
        if((TrackSettingInfo.m_WorkingMode % 10) == WMODE_PERSON)
        {
            int leftShoulderOffset,rightShoulderOffset;
            findShoulderBorder(pu8DiffAbs, IMAGE_WIDTH, IMAGE_HEIGHT,
                m_cx-m_cw/2, m_cy - m_ch/2, m_cx+m_cw/2, m_cy + m_ch/2,
                leftShoulderOffset, rightShoulderOffset);
            float fOffsetRatio = (float)MAX(leftShoulderOffset , rightShoulderOffset)/(m_cw);

            if(fOffsetRatio > 1.0f)
			{
				goto RESET_FAIL;
			}
        }
    }

    printf("reset success bSemiAuto=%d\n", bSemiAuto);
    atRect[1].l32Left = m_cx - m_cw / 2;
    atRect[1].l32Top = m_cy - m_ch / 2;
    atRect[1].l32Width = m_cw;
    atRect[1].l32Height = m_ch;

    //semiauto点到静止区域，直接使用传入的位置
    if(bSemiAuto == TRUE)
    {
        l32 l32CenterShiftX = m_cx - l32InputCenterX;
        l32 l32CenterShiftY = m_cy - l32InputCenterY;
        l32CenterShiftX = ABS(l32CenterShiftX);
        l32CenterShiftY = ABS(l32CenterShiftY);
        if(l32CenterShiftX > m_cw/2 || l32CenterShiftY > m_ch/2)
        {
            //Meanshift之后目标框已经不包含用户所给的点，说明已经偏移过远，此时恢复原始用户输入位置
            m_cw = atRect[0].l32Width;
            m_ch = atRect[0].l32Height;
            m_cx = atRect[0].l32Left + m_cw / 2;
            m_cy = atRect[0].l32Top + m_ch / 2;
        }
    }

    m_f32DiffAvg = 0;
    pu8Src = pu8MaskGmm + (((int)(m_cy - m_ch / 2))*IMAGE_WIDTH);
    u8 *pu8SrcRoi = pu8RoiMask + (((int)(m_cy - m_ch / 2))*IMAGE_WIDTH);
    int roiCnt = 0;
    int i0 = m_cy - m_ch / 2;
    int i1 = m_cy + m_ch / 2;
    int j0 = m_cx - m_cw / 2;
    int j1 = m_cx + m_cw / 2;
    for(int i = i0; i < i1; i++)
    {
        for(int j = j0; j < j1; j++)
        {
            m_f32DiffAvg += pu8Src[j];
            roiCnt += pu8SrcRoi[j] > 0 ? 0 : 1;
        }
        pu8SrcRoi += IMAGE_WIDTH;
        pu8Src+=IMAGE_WIDTH;
    }
    m_f32DiffAvg /= (i1 - i0)*(j1 - j0);


    float f32Wnorm = VanishInfo.GetNormW(m_cw, m_cy + (m_ch/2));
    float f32Hnorm = VanishInfo.GetNormH(m_ch, m_cy + (m_ch/2));

#define BORDER_MARGIN   5
    if(
        (m_cy + (m_ch/2) < IMAGE_HEIGHT - fhuman_h || bSemiAuto == TRUE)//半自动跟踪目标允许跟踪边界位置的目标
        && (m_cx - (m_cw/2) > fhuman_w || bSemiAuto == TRUE)
        && (m_cx + (m_cw/2) < IMAGE_WIDTH - fhuman_w || bSemiAuto == TRUE)
        //&& m_cy - (m_ch/2) > VanishInfo.f32Y0
        && f32Wnorm > 0 && f32Hnorm > 0
        && m_cw > 6 && m_ch > 12
        && (m_f32DiffAvg > 0 || bSemiAuto == TRUE)//半自动跟踪目标允许跟踪静止区域
        /*&& roiCnt ==0*/)
    {       
		//红色框代表meanshift之前目标位置，白色框代表之后目标位置，由于跟踪人时会修改外部传入的尺寸
		//所以看到目标框并不完全和栅栏大小一致
        ShowRgbRec2("meanshift", pu8DiffAbs, IMAGE_WIDTH, IMAGE_HEIGHT, 1, 2, atRect);
        m_bIsSemiAuto = bSemiAuto == TRUE ? 1 : 0;

        //根据是否半自动目标修正目标类型
        if(m_bIsSemiAuto)
        {
            m_iTargetType |= 0x08;
            if(m_f32DiffAvg < 120)
                m_l32MaxTrackTimeMS = g_GblobalSettings.TrackSettingInfo.m_iMaxFollowDuration;
            else
                m_l32MaxTrackTimeMS = 0x7fffffff;
        }
        else
        {
            m_iTargetType &= ~0x08;
            m_l32MaxTrackTimeMS = g_GblobalSettings.TrackSettingInfo.m_iMaxFollowDuration;
        }

        m_l32Life = 0;
        m_bIsMissing = FALSE;

        //检查新目标与前一个跟踪目标是否是同一个目标，不同目标时ID增加
        int tcw = MIN((m_cx + m_cw/2), (fold_cx + fold_cw/2)) - MAX((m_cx - m_cw/2), (fold_cx - fold_cw/2));
        int tch = MIN((m_cy + m_ch/2), (fold_cy + fold_ch/2)) - MAX((m_cy - m_ch/2), (fold_cy - fold_ch/2));
        if((tcw < 0 || tch < 0))
		{
            m_targetid++;

			//计算目标归一化尺寸
			m_Wnorm = f32Wnorm;
			m_Hnorm = f32Hnorm;
		}
		else
		{
			// 同一个目标，如果目标尺寸变化太大，仍然保持原有目标宽度和高度
			if(m_cw < fold_cw * 1.1)
			{
				m_Wnorm = f32Wnorm;
			}
			else
			{
				m_cw = fhuman_w;
			}

			if(m_ch < fold_ch * 1.1)
			{
				m_Hnorm = f32Hnorm;
			}
			else
			{
				m_ch = fhuman_h;
			}
		}

        KalmanFilter8x3::reset(m_cx, m_cy, 0, 0, 0);

        m_QLLSx.reset();
        m_QLLSy.reset();
        m_QLLSx.newData(m_cx);
        m_QLLSy.newData(m_cy);

        // 记录第一次的中心点位置
        m_cx_init = m_cx;
        m_cy_init = m_cy;

		m_bIsTrkHead = FALSE;

		//记录当前目标的中心点位置
		m_last_cy = m_cy;

        m_bIsZoombie = FALSE;
        m_nZoombieCounter = 0;
        m_l32DiffAbnormalCounter = 0;
		m_l32DiffZeroCounter = 0;

        CheckIsZoombie();

        return TRUE;
    }
    else
    {
        printf("reset confirm fail: bottomy = %f, topy = %f, left %f right %f cw = %f, ch = %f human w %f h %f\n",
            m_cy + (m_ch/2), m_cy - (m_ch/2), m_cx - (m_cw/2), m_cx + (m_cw/2),  m_cw, m_ch, fhuman_w, fhuman_h);
    }

    // 		if(checkTargetValid(pu8MaskGmm))
    // 		{
    // 			m_l32Life = 0;
    // 			m_bIsMissing = FALSE;
    // 			m_targetid++;
    // 			KalmanFilter8x3::reset(m_cx, m_cy, 0, 0, 0);
    // 			// 记录第一次的中心点位置
    // 			m_cx_init = m_cx;
    // 			m_cy_init = m_cy;
    //
    // 			return TRUE;
    // 		}
RESET_FAIL:
    if(!m_bIsMissing)
    {
        //恢复之前跟踪的目标
        m_cx = fold_cx;
        m_cy = fold_cy;
        m_cw = fold_cw;
        m_ch = fold_ch;
        m_f32DiffAvg = fold_DiffAvg;
    }
    else
    {
        m_l32Life = -1;
        m_bIsMissing = TRUE;
    }

    return FALSE;
}
BOOL TrackingTarget::CheckIsZoombie()
{
    //如果跟踪框尺寸变得太大,说明已经进入一种不可靠的跟踪状态(Zoombie状态)
    GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;
    if(!m_bIsZoombie)
    {
        f32 f32NormW = VanishInfo.GetNormW(m_cw, m_cy + (m_ch/2));
        if(f32NormW > VanishInfo.f32HumanNormW * 4)
        {
            m_nZoombieCounter ++;
        }
        else
        {
            m_nZoombieCounter = 0;
        }

        if(m_nZoombieCounter > 25)
        {
            m_bIsZoombie = TRUE;
        }
    }
    return m_bIsZoombie;
}
BOOL TrackingTarget::checkTargetValid(u8 *pu8Mask)
{
    f32 f32DiffAvg = 0;
    u8 *pu8OrgMask = pu8Mask;
    GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;

    int y0 = ((int)(m_cy - m_ch / 2));
    int y1 = ((int)(m_cy + m_ch / 2));
    int x0 = ((int)(m_cx - m_cw / 2));
    int x1 = ((int)(m_cx + m_cw / 2));
    if(y0 < 0) y0 = 0;
    if(y1 > IMAGE_HEIGHT) y1 = IMAGE_HEIGHT;
    if(x0 < 0) x0 = 0;
    if(x1 > IMAGE_WIDTH) x1 = IMAGE_WIDTH;

    pu8Mask += y0 * IMAGE_WIDTH;
    for(int i = y0; i < y1; i++)
    {
        for(int j = x0; j < x1; j++)
        {
            f32DiffAvg += pu8Mask[j];
        }
        pu8Mask += IMAGE_WIDTH;
    }
    f32DiffAvg /= (x1 - x0)*(y1 - y0);

    // 		if(m_frameid % 25 == 0)
    // 		{
    // 			UmpPrintf("f32DiffAvg = %f, m_f32DiffAvg = %f\n", f32DiffAvg, m_f32DiffAvg);
    // 			printf("f32DiffAvg = %f, m_f32DiffAvg = %f\n", f32DiffAvg, m_f32DiffAvg);
    // 		}

    m_frameid++;

    // 与目标第一次被选中时候的中心点的距离
    int HaltFlag = 0;
    if((m_l32Life >= 0) && (m_l32Life % 25) == 24)
    {
        float f32dist = sqrt((m_cy - m_cy_init) *(m_cy - m_cy_init) + (m_cx - m_cx_init)*(m_cx - m_cx_init));
        //printf("id = %d, f32dist=%f\n", m_targetid, f32dist);
        if(f32dist < 1)
        {
            HaltFlag = 1;
        }
        m_cx_init = m_cx;
        m_cy_init = m_cy;
    }

    if(f32DiffAvg * 3 < m_f32DiffAvg)
    {
        m_l32DiffAbnormalCounter++;
    }
    else
    {
        m_l32DiffAbnormalCounter = 0;
    }

	if(f32DiffAvg < 0.1)
	{
		m_l32DiffZeroCounter++;
	}
	else
	{
		m_l32DiffZeroCounter = 0;
	}

    //宽高上限
    if(m_cy + (m_ch/2) > IMAGE_HEIGHT
        //|| m_cy - (m_ch/2) < VanishInfo.f32Y0
        || m_cw < 6 || m_ch < 12
        || m_l32DiffAbnormalCounter > 12  // 容易导致目标被删除
        || (m_l32DiffZeroCounter > 3 && (!m_bIsSemiAuto)) //半自动跟踪目标不会因为前景点为全零而消失
        //|| (f32DiffAvg == 0 && (!m_bIsSemiAuto)) //半自动跟踪目标不会因为前景点为全零而消失
        || (HaltFlag == 1 && (!m_bIsSemiAuto)))  //半自动跟踪目标不会因为静止而消失
    {
        printf("checkvalid :id = %d, m_cy = %f cw = %f, ch = %f, DifAvg=(%f/%f)  HaltFlag=%d m_l32DiffAbnormalCounter %d SemiAuto %d\n",
            m_targetid, m_cy, m_cw, m_ch,
            f32DiffAvg, m_f32DiffAvg,
            HaltFlag, m_l32DiffAbnormalCounter, m_bIsSemiAuto
            );
        return FALSE;
    }
    return TRUE;
}

static BOOL IsSameRect(const POINT& tTL1, const POINT& tBR1, const POINT& tTL2, const POINT& tBR2)
{
	l32 l32Area1;
	l32 l32Area2;

	l32 l32Left = MAX(tTL1.x, tTL2.x);
	l32 l32Right = MIN(tBR1.x, tBR2.x);
	l32 l32Top = MAX(tTL1.y, tTL2.y);
	l32 l32Bottom = MIN(tBR1.y, tBR2.y);
	if(l32Right <= l32Left || l32Bottom <= l32Top)
		return FALSE;

	l32 l32CommonArea = (l32Right - l32Left + 1) * (l32Bottom - l32Top + 1);
	l32Area1 = (tBR1.x - tTL1.x + 1) * (tBR1.y - tTL1.y + 1);
	l32Area2 = (tBR2.x - tTL2.x + 1) * (tBR2.y - tTL2.y + 1);
	if(l32CommonArea > 6 * l32Area1 / 10 && l32CommonArea > 6 * l32Area2 / 10)
		return TRUE;

	return FALSE;
}

static BOOL IsSameRectTrk(const POINT& tTL1, const POINT& tBR1, const POINT& tTL2, const POINT& tBR2, l32 l32Mass)
{
	l32 l32Area1;
	l32 l32Area2;

	l32 l32Left = MAX(tTL1.x, tTL2.x);
	l32 l32Right = MIN(tBR1.x, tBR2.x);
	l32 l32Top = MAX(tTL1.y, tTL2.y);
	l32 l32Bottom = MIN(tBR1.y, tBR2.y);
	if(l32Right <= l32Left || l32Bottom <= l32Top)
		return FALSE;

	l32 l32CommonArea = (l32Right - l32Left + 1) * (l32Bottom - l32Top + 1);
	l32Area1 = (tBR1.x - tTL1.x + 1) * (tBR1.y - tTL1.y + 1);
	l32Area2 = (tBR2.x - tTL2.x + 1) * (tBR2.y - tTL2.y + 1);

	if((l32CommonArea > 6 * l32Area1 / 10 && l32CommonArea > 6 * l32Area2 / 10)) // || (l32Mass > l32Area1 * 65 / 100))
		return TRUE;

	return FALSE;
}


BOOL TrackingTarget::FindRectByTrack(f32 &f32CurX, f32 &f32CurY,
									 f32 &l32WndW, f32 &l32WndH,
									 f32 f32PreX, f32 f32PreY)
{
	GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;
	GTargetInfo &TargetInfo = g_GblobalSettings.TargetInfo;

	BOOL bFind = FALSE;

	l32 l32Idx, l32IntersectIdx, l32IntersectIdx1, l32IntersectNum;
	float f32Trackw;
	l32 l32Left, l32Right, l32Top, l32Bottom;
	l32 l32CommonArea, l32MaxCommonArea = 0;
	
	POINT tTL, tBR;
	
	tTL.x = f32PreX - l32WndW / 2;
	tTL.y = f32PreY - l32WndH / 2;
	tBR.x = f32PreX + l32WndW / 2;
	tBR.y = f32PreY + l32WndH / 2;


	l32IntersectNum = 0;

	for(l32Idx = 0; l32Idx < TargetInfo.l32TrackNum; l32Idx++)
	{
		//
		l32Left = MAX(tTL.x, TargetInfo.atPointL[l32Idx].x);
		l32Right = MIN(tBR.x, TargetInfo.atPointR[l32Idx].x);
		l32Top = MAX(tTL.y, TargetInfo.atPointL[l32Idx].y);
		l32Bottom = MIN(tBR.y, TargetInfo.atPointR[l32Idx].y);

		//非相交
		if(l32Right <= l32Left || l32Bottom <= l32Top)
		{
			continue;
		}

		l32CommonArea = (l32Right - l32Left + 1) * (l32Bottom - l32Top + 1);

		//目标上半身宽度在人的宽度[0.5,1.5]之间，
		f32Trackw = VanishInfo.GetW(VanishInfo.f32HumanNormW, TargetInfo.atPointR[l32Idx].y);

		l32IntersectNum++;
		l32IntersectIdx1 = l32Idx;

		//相加面积最大，且两个目标尺寸相差不大，
		if(l32CommonArea > l32MaxCommonArea && (TargetInfo.al32UperBodyW[l32Idx] < f32Trackw * 2.0))
		{
			l32MaxCommonArea = l32CommonArea;
			l32IntersectIdx = l32Idx;
		}
	}

	//直接替换条件：当前目标与跟踪目标只有一个有相交，该目标一对一匹配，两个目标尺寸相差不大
	if((1 == l32IntersectNum && TargetInfo.al321v1Cnt[l32IntersectIdx1] > 2))
	{
		f32CurX = (TargetInfo.atPointL[l32IntersectIdx1].x + TargetInfo.atPointR[l32IntersectIdx1].x) / 2;
		f32CurY = (TargetInfo.atPointL[l32IntersectIdx1].y + TargetInfo.atPointR[l32IntersectIdx1].y) / 2;

		bFind = TRUE;
	}
	else if(l32IntersectNum > 1 && l32MaxCommonArea > 0)
	{
		if(TargetInfo.al321v1Cnt[l32IntersectIdx] > 2)
		{
			f32CurX = (TargetInfo.atPointL[l32IntersectIdx].x + TargetInfo.atPointR[l32IntersectIdx].x) / 2;
			f32CurY = (TargetInfo.atPointL[l32IntersectIdx].y + TargetInfo.atPointR[l32IntersectIdx].y) / 2;

			bFind = TRUE;
		}
	}

	return bFind;

}

int x,y,w,h;
void TrackingTarget::tracking_byKLT(u8 *pu8Pre, u8 *pu8Cur, l32 l32Width, l32 l32Height, u8 *pu8WeightMap, u8 *pu8Mask, AntBlob **blobs, int nBlobs)
{
    float f_error, f_nx, f_ny, f_mass;
    float f_nw, f_nh;

    GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;

    if(m_bIsMissing) return;
    //	if(m_l32Life < 0) return;

    //Kalman预测新位置
    f_nx = m_cx;
    f_ny = m_cy;
    KalmanFilter8x3::predictState(f_nx, f_ny, f_mass);    

    //以人的尺度为跟踪框限制,因此即使是出现人群,也仅仅以其中接近画面中心位置的一个人为主要跟踪目标
    f_nw = VanishInfo.GetW((m_Wnorm + VanishInfo.f32HumanNormW) * 0.5f, f_ny + (m_ch/2));
	f_nh = VanishInfo.GetH(VanishInfo.f32HumanNormH, f_ny + (m_ch/2));

    //光流法内部资源限制,不能跟踪太大的Patch(在704x396图像上很少有目标需要这么大的跟踪区域)
    if(f_nw > 64) f_nw = 64;
    if(f_nh > 64) f_nh = 64;
	float f_nx_rs = f_nx;
	float f_ny_rs = f_ny;
    KLTTrackingPointWithWeight(pu8Pre, pu8Cur, pu8WeightMap, l32Width, l32Height,
        m_cx, m_cy, f_nx, f_ny,
        f_nw, f_nh, &f_error);
    // 光流法之后，依据findblob的结果对跟踪框的大小和位置进行修正

	x = m_cx;
	y = m_cy;
	w = m_cw;
	h = m_ch;
#if 1 // 背景建模与光流结合
	f32 cx_gmm = f_nx;
	f32 cy_gmm = f_ny;
	f32 cw_gmm = m_cw;
	f32 ch_gmm = m_ch;

	BOOL bRet = FindRectByTrack(cx_gmm, cy_gmm, cw_gmm, ch_gmm, m_cx, m_cy);
	if(bRet)
	{
		x = f_nx;
		y = f_ny;
		f_ny = cy_gmm;
		f_nx = cx_gmm;
	}
	else
	{
		x = -1;
	}


	KalmanFilter8x3::updateState(f_nx, f_ny, f_mass);
	m_cx = f_nx;
	m_cy = f_ny;

	m_cw = VanishInfo.GetW(m_Wnorm, m_cy + (m_ch/2));//m_Wnorm * (m_cy + (m_ch/2) - VanishInfo.f32Y0);
	m_ch = VanishInfo.GetH(m_Hnorm, m_cy + (m_ch/2));//m_Hnorm * (m_cy + (m_ch/2) - VanishInfo.f32Y0);

	w = m_cw;
	h = m_ch;

#endif
    //printf("x4 %f,%f   %fx%f\n",f_nx,f_ny,m_cw,m_ch);
    //if(m_ch > (g_MaxTargetHeight * (m_cy + (m_ch/2) - g_fVanishY0)/(g_TargetY - g_fVanishY0)))
    //	printf(">>>>>>>>>>>>>>>>>>>>>>>>>m_ch > g_MaxTargetHeight\n");

    //记录位置到预测器中
    m_QLLSx.newData(m_cx);
    m_QLLSy.newData(m_cy);

	m_last_cy = m_cy;

    CheckIsZoombie();

    if(!checkTargetValid(pu8Mask))
    {
        m_bIsMissing = TRUE;
        m_l32Life = -1;
        m_f32DiffAvg = 0.0f;
    }
    else
    {
        m_bIsMissing = FALSE;
        m_l32Life++;
    }
}

BOOL TrackingTarget::GetCurTargetRect(TRect * ptRect)
{
    if(!m_bIsMissing && !m_bIsZoombie)
    {
        ptRect->l32Left = m_cx - m_cw/2;
        ptRect->l32Top = m_cy - m_ch/2;
        ptRect->l32Width = m_cw;
        ptRect->l32Height = m_ch;
        return TRUE;
    }
    return FALSE;
}
#include <stdio.h>
#include <string.h>
#include "ai_defs.h"
#include <math.h>

typedef float f32;
#define OPENCV_FLT_EPSILON     1.192092896e-07F        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
#define VECTHRESH    0.1  //运动矢量的阈值
#define MAX_RCID_CNT_PER_GROUP  100
#define MAX_WNDSIZE    70
#define MAXITERATION   30

#ifdef WIN32
#include "intrinsic.h"
#else
#endif

/*=============================================================================
函 数 名: KLTGetSubImage
功    能: 提取目标子图像
算法实现: 
全局变量: 无
参    数:   pu8Src		  当前帧灰度图像[in]
l32Width      当前帧宽度[in]
l32Height     当前帧高度[in]
pu8Dst        提取目标图像[out]
l32W          目标图像所在矩形框的宽度[in]
l32H          目标图像所在矩形框的高度[in]
f32X          目标图像特征l32X方向上的运动矢量[in]
f32Y          目标图像特征点l32Y方向上的运动矢量[in]
返 回 值: 无
=============================================================================*/
static void KLTGetSubImage(u8 *pu8Src, l32 l32Width, l32 l32Height,
                           u8 *pu8Dst, l32 l32W, l32 l32H,
                           f32 f32X, f32 f32Y)
{
    l32 l32X, l32Y , l32Sx, l32Sy, l32Cx, l32Cy, l32Sx2;
    l32 l32Val0, l32Val1, l32Wx, l32Wy;
    u8 *pu8Row0;
    u8 *pu8Row1;

    //需要提取图像窗口中心坐标的整数部分，小数部分
    l32Cx = (l32)f32X;
    l32Cy = (l32)f32Y;
    l32Wx = (l32)((f32X - l32Cx) * 16);
    l32Wy = (l32)((f32Y - l32Cy) * 16);

    //左上角点坐标
    l32Cx -= (l32W >> 1);
    l32Cy -= (l32H >> 1);
    for(l32Y = 0; l32Y < l32H; l32Y++)
    {
        l32Sy = l32Y + l32Cy;
        if(l32Sy < 0)
        {
            l32Sy=0;
        }
        if(l32Sy > l32Height - 1)
        {
            l32Sy = l32Height - 1;
        }

        pu8Row0 = pu8Src + l32Sy * l32Width;
        if(l32Sy == l32Height - 1)
            pu8Row1 = pu8Row0;
        else
            pu8Row1 = pu8Row0 + l32Width;
        for(l32X = 0; l32X < l32W; l32X++)
        {
            //源图对应整像素点左上角坐标(l32Sx,l32Sy)
            l32Sx = l32X + l32Cx;
            if(l32Sx < 0)
            {
                l32Sx = 0;
            }
            l32Sx2 = l32Sx + 1;
            if(l32Sx > l32Width - 1)
            {
                l32Sx2 = l32Sx = l32Width - 1;
            }

            //插值
            l32Val0 = pu8Row0[l32Sx] * (16 - l32Wx) + pu8Row0[l32Sx2] * (l32Wx);
            l32Val1 = pu8Row1[l32Sx] * (16 - l32Wx) + pu8Row1[l32Sx2] * (l32Wx);
            pu8Dst[l32X] = (u8)((l32Val0 * (16 - l32Wy) + l32Val1 * l32Wy) / 256);
        }
        pu8Dst += l32W;
    }
}

/*=============================================================================
函数名：KLTTrackingPoint
功   能：KL迭代的特征点跟踪
算法实现：
sum(f32Ix*f32Ix)    sum(f32Ix*f32Iy)
G= [                                     ]  其中sum是在窗口最大为21x21的临域求和
sum(f32Ix*f32Iy)    sum(f32Iy*f32Iy)     
sum([A(x,y) - B(x,y)] * f32Ix)
b = [                             ]
sum([A(x,y) - B(x,y)] * f32Iy)
v = revers(G) * b
全局变量：无
参   数：pu8Pre        前一帧灰度图像[in]
pu8Cur        当前帧灰度图像[in]
l32Width      图像宽度[in]
l32Height     图像高度[in]
f32PrevX      目标在前一帧中x方向上的运动矢量[in]
f32PrevY      目标在前一帧中y方向上的运动矢量[in]
f32CurX       目标在当前帧中x方向上的运动矢量[out]
f32CurY       目标在当前帧中y方向上的运动矢量[out]
l32WndW       目标所在矩形框的宽度[in]
l32WndH       目标所在矩形框的高度[in]

返回值：是否跟踪成功
=============================================================================*/
l32 KLTTrackingPoint(u8 *pu8Pre, u8 *pu8Cur, l32 l32Width, l32 l32Height,
                     f32 f32PrevX, f32 f32PrevY, f32 &f32CurX, f32 &f32CurY,
                     l32 l32WndW, l32 l32WndH, f32 *pf32Error)
{
    l32 l32I, l32X, l32Y, l32Count;
    l32 l32Cx, l32Cy;
    f32 f32Ix, f32Iy;
    l32 l32Ix, l32Iy, l32II;
    f32 f32Gxx, f32Gxy, f32Gyy, f32Det, f32Ex, f32Ey, f32D;
    f32 f32Vx, f32Vy;

    u8 au8A[80 * 80] = {0};
    u8 au8B[MAX_WNDSIZE * MAX_WNDSIZE] = {0};
    s16 as16D[MAX_WNDSIZE * MAX_WNDSIZE * 2] = {0};
    u8 *pu8A;
    u8 *pu8B;
    l32 l32Res = 1;

    //提取前一帧中某个矩形框内的目标图像
    //因为后续要对目标图像求梯度，所以矩形框左右，上下各向外扩1个像素
    KLTGetSubImage(pu8Pre, l32Width, l32Height,
        au8A, l32WndW + 2,  l32WndH + 2, f32PrevX, f32PrevY);

    f32Gxx = f32Gxy = f32Gyy = 0.0f;

    //对目标图像A求梯度
    l32I = 0;
    for(l32Y = 1; l32Y <= l32WndH; l32Y++)
    {
        pu8A = au8A + (l32WndW + 2) * l32Y;
        for(l32X=1; l32X<= l32WndW; l32X++)
        {
            l32Ix = ((l32)pu8A[l32X + 1]         - (l32)pu8A[l32X - 1])/2;
            l32Iy = ((l32)pu8A[l32X + (l32WndW + 2)] - (l32)pu8A[l32X-(l32WndW + 2)])/2;
            f32Gxx += l32Ix * l32Ix;
            f32Gxy += l32Ix * l32Iy;
            f32Gyy += l32Iy * l32Iy;
            as16D[l32I * 2] = l32Ix;   //保存当前像素点在x方向上的梯度
            as16D[l32I * 2 + 1] = l32Iy; //保存当前像素点在y方向上的梯度
            l32I++;
        }
    }

    f32Det = f32Gxx * f32Gyy - f32Gxy * f32Gxy;
    if (f32Det < OPENCV_FLT_EPSILON)
    {
        l32Res = 0;
        return l32Res;
    }//方程的解不稳定，认为跟踪失败

    f32Det = 1/f32Det;
    f32Gxx *= f32Det;
    f32Gyy *= f32Det;
    f32Gxy *= f32Det;

    for(l32Count = 0; l32Count < MAXITERATION; l32Count++)//迭代多次，得到目标精确的运动矢量
    {
        l32I = 0;
        f32Ex = f32Ey = 0.0f;
        //更新目标图B(每次都要重新插值)
        KLTGetSubImage(pu8Cur, l32Width, l32Height,
            au8B, l32WndW, l32WndH, f32CurX, f32CurY);
        for(l32Y = 1; l32Y <= l32WndH; l32Y++)
        {
            pu8A = au8A + (l32WndW + 2) * l32Y;
            pu8B = au8B + l32WndW * (l32Y - 1);
            for(l32X = 1; l32X <= l32WndW; l32X++)
            {
                l32II = pu8A[l32X] - pu8B[l32X - 1];
                f32Ex += l32II * as16D[l32I * 2];
                f32Ey += l32II * as16D[l32I * 2 + 1];
                l32I++;
            }
        }
        f32Vx = (f32Gyy * f32Ex - f32Gxy * f32Ey);
        f32Vy = (f32Gxx * f32Ey - f32Gxy * f32Ex);
        f32CurX += f32Vx;
        f32CurY += f32Vy;
        if(fabs(f32Vx) < VECTHRESH && fabs(f32Vy) < VECTHRESH)
            break;
    }
    //printf("KLT iterate %d times\n", l32Count);

    //计算误差
    if(pf32Error)
    {
        *pf32Error = 0;
        for(l32Y = 1; l32Y <= l32WndH; l32Y++)
        {
            pu8A = au8A + (l32WndW + 2) * l32Y;
            pu8B = au8B + l32WndW * (l32Y - 1);
            for(l32X = 1; l32X <= l32WndW; l32X++)
            {
                f32D = pu8A[l32X] - pu8B[l32X - 1];
                *pf32Error += fabs(f32D);
            }
        }
    }

    return l32Res;
}

 l32 KLTTrackingPointWithWeight2(u8 *pu8Pre, u8 *pu8Cur, u8 *pu8PreWeight, l32 l32Width, l32 l32Height,
	 f32 f32PrevX, f32 f32PrevY, f32 &f32CurX, f32 &f32CurY,
	 l32 l32WndW, l32 l32WndH, f32 *pf32Error)
 {
	 l32 l32I, l32X, l32Y, l32Count;
	 l32 l32Cx, l32Cy;
	 l32 l32Ix, l32Iy, l32II, l32W;
	 f32 f32Ix, f32Iy;
	 f32 f32Gxx, f32Gxy, f32Gyy, f32Det, f32Ex, f32Ey, f32D;
	 f32 f32Vx, f32Vy;

	 u8 au8A[80 * 80] = {0};
	 u8 au8W[80 * 80] = {0};
	 u8 au8B[MAX_WNDSIZE * MAX_WNDSIZE] = {0};
	 l32 al32D[MAX_WNDSIZE * MAX_WNDSIZE * 2] = {0};
	 u8 *pu8A;
	 u8 *pu8W;
	 u8 *pu8B;
	 l32 l32Res = 1;

	 //提取前一帧中某个矩形框内的目标图像
	 //因为后续要对目标图像求梯度，所以矩形框左右，上下各向外扩1个像素
	 KLTGetSubImage(pu8Pre, l32Width, l32Height,
		 au8A, l32WndW + 2,  l32WndH + 2, f32PrevX, f32PrevY);

	 KLTGetSubImage(pu8PreWeight, l32Width, l32Height,
		 au8W, l32WndW + 2,  l32WndH + 2, f32PrevX, f32PrevY);

	 f32Gxx = f32Gxy = f32Gyy = 0.0f;

	 //对目标图像A求梯度
	 l32I = 0;
	 for(l32Y = 1; l32Y <= l32WndH; l32Y++)
	 {
		 pu8A = au8A + (l32WndW + 2) * l32Y;
		 pu8W = au8W + (l32WndW + 2) * l32Y;
		 for(l32X=1; l32X<= l32WndW; l32X++)
		 {
			 l32Ix = ((l32)pu8A[l32X + 1]         - (l32)pu8A[l32X - 1])/2;
			 l32Iy = ((l32)pu8A[l32X + (l32WndW + 2)] - (l32)pu8A[l32X-(l32WndW + 2)])/2;
			 l32W = pu8W[l32X];//权重

			 f32Gxx += l32W * l32Ix * l32Ix;
			 f32Gxy += l32W * l32Ix * l32Iy;
			 f32Gyy += l32W * l32Iy * l32Iy;
			 al32D[l32I * 2] = l32W * l32Ix;   //保存当前像素点在x方向上的梯度
			 al32D[l32I * 2 + 1] = l32W * l32Iy; //保存当前像素点在y方向上的梯度
			 l32I++;
		 }
	 }

	 f32Det = f32Gxx * f32Gyy - f32Gxy * f32Gxy;
	 if (f32Det < OPENCV_FLT_EPSILON)
	 {
		 l32Res = 0;
		 return l32Res;
	 }//方程的解不稳定，认为跟踪失败

	 for(l32Count = 0; l32Count < MAXITERATION; l32Count++)//迭代多次，得到目标精确的运动矢量
	 {
		 l32I = 0;
		 f32Ex = f32Ey = 0.0f;
		 //更新目标图B(每次都要重新插值)
		 KLTGetSubImage(pu8Cur, 704, 396,
			 au8B, l32WndW, l32WndH, f32CurX, f32CurY);
		 for(l32Y = 1; l32Y <= l32WndH; l32Y++)
		 {
			 pu8A = au8A + (l32WndW + 2) * l32Y;
			 pu8B = au8B + l32WndW * (l32Y - 1);
			 for(l32X = 1; l32X <= l32WndW; l32X++)
			 {
				 l32II = pu8A[l32X] - pu8B[l32X - 1];
				 f32Ex += l32II * al32D[l32I * 2];
				 f32Ey += l32II * al32D[l32I * 2 + 1];
				 l32I++;
			 }
		 }
		 f32Vx = (f32Gyy * f32Ex - f32Gxy * f32Ey) / f32Det;
		 f32Vy = (f32Gxx * f32Ey - f32Gxy * f32Ex) / f32Det;
		 f32CurX += f32Vx;
		 f32CurY += f32Vy;
		 if(fabs(f32Vx) < VECTHRESH && fabs(f32Vy) < VECTHRESH)
			 break;
	 }
	 //printf("KLT iterate %d times\n", l32Count);

	 //计算误差
	 if(pf32Error)
	 {
		 *pf32Error = 0;
		 for(l32Y = 1; l32Y <= l32WndH; l32Y++)
		 {
			 pu8A = au8A + (l32WndW + 2) * l32Y;
			 pu8B = au8B + l32WndW * (l32Y - 1);
			 for(l32X = 1; l32X <= l32WndW; l32X++)
			 {
				 f32D = pu8A[l32X] - pu8B[l32X - 1];
				 *pf32Error += fabs(f32D);
			 }
		 }
	 }

	 return l32Res;
}


l32 KLTTrackingPointWithWeight(u8 *pu8Pre, u8 *pu8Cur, u8 *pu8PreWeight, l32 l32Width, l32 l32Height,
                               f32 f32PrevX, f32 f32PrevY, f32 &f32CurX, f32 &f32CurY,
                               l32 l32WndW, l32 l32WndH, f32 *pf32Error)
{
    l32 l32I, l32X, l32Y, l32Count;
    l32 l32Cx, l32Cy;
    l32 l32Ix, l32Iy, l32II, l32W;
    f32 f32Ix, f32Iy;
    f32 f32Gxx, f32Gxy, f32Gyy, f32Det, f32Ex, f32Ey, f32D;
    f32 f32Vx, f32Vy;

    u8 au8A[80 * 80] = {0};
    u8 au8W[80 * 80] = {0};
    u8 au8B[MAX_WNDSIZE * MAX_WNDSIZE] = {0};
    l32 al32D[MAX_WNDSIZE * MAX_WNDSIZE * 2] = {0};
    u8 *pu8A;
    u8 *pu8W;
    u8 *pu8B;
    l32 l32Res = 1;

    //提取前一帧中某个矩形框内的目标图像
    //因为后续要对目标图像求梯度，所以矩形框左右，上下各向外扩1个像素
    KLTGetSubImage(pu8Pre, l32Width, l32Height,
        au8A, l32WndW + 2,  l32WndH + 2, f32PrevX, f32PrevY);

    KLTGetSubImage(pu8PreWeight, l32Width, l32Height,
        au8W, l32WndW + 2,  l32WndH + 2, f32PrevX, f32PrevY);

    f32Gxx = f32Gxy = f32Gyy = 0.0f;

    //对目标图像A求梯度
    l32I = 0;
    for(l32Y = 1; l32Y <= l32WndH; l32Y++)
    {
        pu8A = au8A + (l32WndW + 2) * l32Y;
        pu8W = au8W + (l32WndW + 2) * l32Y;
        for(l32X=1; l32X<= l32WndW; l32X++)
        {
            l32Ix = ((l32)pu8A[l32X + 1]         - (l32)pu8A[l32X - 1])/2;
            l32Iy = ((l32)pu8A[l32X + (l32WndW + 2)] - (l32)pu8A[l32X-(l32WndW + 2)])/2;
            l32W = pu8W[l32X];//权重

            f32Gxx += l32W * l32Ix * l32Ix;
            f32Gxy += l32W * l32Ix * l32Iy;
            f32Gyy += l32W * l32Iy * l32Iy;
            al32D[l32I * 2] = l32W * l32Ix;   //保存当前像素点在x方向上的梯度
            al32D[l32I * 2 + 1] = l32W * l32Iy; //保存当前像素点在y方向上的梯度
            l32I++;
        }
    }

    f32Det = f32Gxx * f32Gyy - f32Gxy * f32Gxy;
    if (f32Det < OPENCV_FLT_EPSILON)
    {
        l32Res = 0;
        return l32Res;
    }//方程的解不稳定，认为跟踪失败

    for(l32Count = 0; l32Count < MAXITERATION; l32Count++)//迭代多次，得到目标精确的运动矢量
    {
        l32I = 0;
        f32Ex = f32Ey = 0.0f;
        //更新目标图B(每次都要重新插值)
        KLTGetSubImage(pu8Cur, l32Width, l32Height,
            au8B, l32WndW, l32WndH, f32CurX, f32CurY);
        for(l32Y = 1; l32Y <= l32WndH; l32Y++)
        {
            pu8A = au8A + (l32WndW + 2) * l32Y;
            pu8B = au8B + l32WndW * (l32Y - 1);
            for(l32X = 1; l32X <= l32WndW; l32X++)
            {
                l32II = pu8A[l32X] - pu8B[l32X - 1];
                f32Ex += l32II * al32D[l32I * 2];
                f32Ey += l32II * al32D[l32I * 2 + 1];
                l32I++;
            }
        }
        f32Vx = (f32Gyy * f32Ex - f32Gxy * f32Ey) / f32Det;
        f32Vy = (f32Gxx * f32Ey - f32Gxy * f32Ex) / f32Det;
        f32CurX += f32Vx;
        f32CurY += f32Vy;
        if(fabs(f32Vx) < VECTHRESH && fabs(f32Vy) < VECTHRESH)
            break;
    }
    //printf("KLT iterate %d times\n", l32Count);

    //计算误差
    if(pf32Error)
    {
        *pf32Error = 0;
        for(l32Y = 1; l32Y <= l32WndH; l32Y++)
        {
            pu8A = au8A + (l32WndW + 2) * l32Y;
            pu8B = au8B + l32WndW * (l32Y - 1);
            for(l32X = 1; l32X <= l32WndW; l32X++)
            {
                f32D = pu8A[l32X] - pu8B[l32X - 1];
                *pf32Error += fabs(f32D);
            }
        }

		*pf32Error /= l32WndH * l32WndW;
    }

    return l32Res;
}



void SmoothFilter7x7(u8 *pu8Src, l32 l32Width, l32 l32Height, u8 *pu8Dst, u8 *pu8Temp7Line)
{
    l32 l32x, l32y, l32id, l32Sum, l32Index;
    u8 * restrict pu8HResult;
    u8 * restrict pu8LineN3;
    u8 * restrict pu8LineN2;
    u8 * restrict pu8LineN1;
    u8 * restrict pu8Line0;
    u8 * restrict pu8LineP1;
    u8 * restrict pu8LineP2;
    u8 * restrict pu8LineP3;
    u64 u64Data0;
    u64 u64DataP1, u64DataP2, u64DataP3;
    u64 u64DataN1, u64DataN2, u64DataN3;

    u32 u32DataL0, u32DataH0;
    u32 u32DataL1, u32DataH1;
    u32 u32Data0, u32Data1, u32Data2, u32Data3, u32Data4;

    memset(pu8Temp7Line, 0, l32Width*7);
    for(l32y = 0, l32id = 0; l32y < l32Height + 2; l32y++, pu8Src+= l32Width)
    {
        pu8HResult = pu8Temp7Line + (l32id) * l32Width;

        if(l32y < l32Height)
        {
            pu8HResult[0] = pu8HResult[1] = pu8HResult[2] = 0;
            l32x = 3;
            //当前行水平7点求和
            for(; l32x < l32Width - 8; l32x += 8)
            {
                u64DataN3 = _mem8_const(pu8Src + l32x - 3);
                u64DataN2 = _mem8_const(pu8Src + l32x - 2);
                u64DataN1 = _mem8_const(pu8Src + l32x - 1);
                u64Data0  = _mem8_const(pu8Src + l32x);
                u64DataP1 = _mem8_const(pu8Src + l32x + 1);
                u64DataP2 = _mem8_const(pu8Src + l32x + 2);
                u64DataP3 = _mem8_const(pu8Src + l32x + 3);

                u32DataL0 = _avgu4(_avgu4(_loll(u64DataP2), _loll(u64DataN2)), _avgu4(_loll(u64DataP3), _loll(u64DataN3)));
                u32DataL1 = _avgu4(_avgu4(_loll(u64DataP1), _loll(u64DataN1)), _loll(u64Data0));
                u32DataL1 = _avgu4(u32DataL0, u32DataL1);

                u32DataH0 = _avgu4(_avgu4(_hill(u64DataP2), _hill(u64DataN2)), _avgu4(_hill(u64DataP3), _hill(u64DataN3)));
                u32DataH1 = _avgu4(_avgu4(_hill(u64DataP1), _hill(u64DataN1)), _hill(u64Data0));
                u32DataH1 = _avgu4(u32DataH0, u32DataH1);

                _mem8(pu8HResult + l32x) = _itoll(u32DataH1, u32DataL1);
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
            l32Index = l32id >= 3? (l32id-3):(l32id + 7 - 3);
            pu8LineN3 = pu8Temp7Line + (l32Index) * l32Width;
            l32Index = l32id >= 2? (l32id-2):(l32id + 7 - 2);
            pu8LineN2 = pu8Temp7Line + (l32Index) * l32Width;
            l32Index = l32id >= 1? (l32id-1):(l32id + 7 - 1);
            pu8LineN1 = pu8Temp7Line + (l32Index) * l32Width;

            pu8Line0 = pu8HResult;

            l32Index = l32id < (7-1)? (l32id+1):(l32id + 1 - 7);
            pu8LineP1 = pu8Temp7Line + (l32Index) * l32Width;
            l32Index = l32id < (7-2)? (l32id+2):(l32id + 2 - 7);
            pu8LineP2 = pu8Temp7Line + (l32Index) * l32Width;
            l32Index = l32id < (7-3)? (l32id+3):(l32id + 3 - 7);
            pu8LineP3 = pu8Temp7Line + (l32Index) * l32Width;

            l32x = 0;
            for(; l32x < l32Width; l32x += 8)
            {
                u64DataN3 = _amem8(pu8LineN3 + l32x);
                u64DataN2 = _amem8(pu8LineN2 + l32x);
                u64DataN1 = _amem8(pu8LineN1 + l32x);
                u64Data0 = _amem8(pu8Line0 + l32x);
                u64DataP1 = _amem8(pu8LineP1 + l32x);
                u64DataP2 = _amem8(pu8LineP2 + l32x);
                u64DataP3 = _amem8(pu8LineP3 + l32x);

                u32DataL0 = _avgu4(_avgu4(_loll(u64DataP2), _loll(u64DataN2)), _avgu4(_loll(u64DataP3), _loll(u64DataN3)));
                u32DataL1 = _avgu4(_avgu4(_loll(u64DataP1), _loll(u64DataN1)), _loll(u64Data0));
                u32DataL1 = _avgu4(u32DataL0, u32DataL1);

                u32DataH0 = _avgu4(_avgu4(_hill(u64DataP2), _hill(u64DataN2)), _avgu4(_hill(u64DataP3), _hill(u64DataN3)));
                u32DataH1 = _avgu4(_avgu4(_hill(u64DataP1), _hill(u64DataN1)), _hill(u64Data0));
                u32DataH1 = _avgu4(u32DataH0, u32DataH1);

                _amem8(pu8Dst + l32x) = _itoll(u32DataH1, u32DataL1);
            }
            pu8Dst += l32Width;
        }
        l32id ++;
        if(l32id >= 7) l32id = 0;
    }
}



static void getII(u32 * ii, u8 *pu8Src, int stride, int w, int h)
{
    //积分图像
    u8 *pu8base = pu8Src;
    u32 * base = ii;
    u32 *prev_line = ii;
    u32 s;

    //第一行，全零
    for (int x = 0; x < w+1; x++) ii[x] = 0;
    ii += (w + 1);
    //每行第一列也是全零
    ii[0] = 0;
    for (int x = 0; x < w; x++) {
        ii[x+1] = pu8Src[x] + ii[x];
    }
    prev_line = ii;
    ii += (w + 1); pu8Src += stride;

    for (int y = 1; y < h; y++) {
        ii[0] = 0;
        s = 0;
        for (int x = 0; x < w; x++) {
            s += pu8Src[x];
            ii[x + 1] = s + prev_line[x+1];
        }
        prev_line = ii;
        ii += (w + 1); pu8Src += stride;
    }
}

void SubAbsImages(u8 * restrict pu8Pre, u8 * restrict pu8Cur, u8 * restrict pu8Dst,
                  int w, int h)
{
    u64 u64Data0, u64Data1;
    u32 u32Data0, u32Data1;
    int x, y;
    for(y=0; y<h; y++)
    {
        for(x=0; x<w; x+=8)
        {
            u64Data0 = _amem8/*_const*/(pu8Cur + x);
            u64Data1 = _amem8/*_const*/(pu8Pre + x);
            u32Data0 = _subabs4(_loll(u64Data0),_loll(u64Data1));
            u32Data1 = _subabs4(_hill(u64Data0),_hill(u64Data1));
            _amem8(pu8Dst + x) = _itoll(u32Data1, u32Data0);
        }
        pu8Cur += w;
        pu8Pre += w;
        pu8Dst += w;
    }
}



#ifdef WIN32
#define OPENCV_DEBUG 1
#else
#define OPENCV_DEBUG 0
#endif

int GetTRectOverlapArea(TRect * pRc1, TRect * pRc2)
{
    int max_xs = MAX(pRc1->l32Left, pRc2->l32Left);
    int min_xe = MIN(pRc1->l32Left + pRc1->l32Width , pRc2->l32Left + pRc2->l32Width);
    int ov_x = min_xe - max_xs;

    int max_ys = MAX(pRc1->l32Top, pRc2->l32Top);
    int min_ye = MIN(pRc1->l32Top + pRc1->l32Height, pRc2->l32Top + pRc2->l32Height);
    int ov_y = min_ye - max_ys;

    if(ov_x < 0 || ov_y < 0) return 0;

    return ov_x * ov_y;
}

static BOOL IsExistFg(TRect *ptObjRc, u8 *pu8FgMask)
{
    BOOL bIsExitFg = 0;
    u8 *pu8Src = pu8FgMask + ptObjRc->l32Top * IMAGE_WIDTH + ptObjRc->l32Left;
    l32 l32Cnt = 0;
    for(int i = 0; i < ptObjRc->l32Height; i++)
    {
        for(int j = 0; j < ptObjRc->l32Width; j++)
        {
            l32Cnt += pu8Src[j] > 0 ? 1 : 0;
        }
        pu8Src+=IMAGE_WIDTH;
    }
    //	if(l32Cnt < (l32)(ptObjRc->l32Width * (float)ptObjRc->l32Height * 0.05f)) bIsExitFg = 1;
    if(l32Cnt > 0) bIsExitFg = 1;
    return bIsExitFg;
}

void DetectMovingHumanByFENCE(u8 * restrict pu8Diff,   //背景减除之后的图像
                              u8 * restrict pu8FgMask, //Gmm前景图像mask
                              TRect * pCurRect,        //当前跟踪目标
                              TRect * pNewRect)
{
    int w, h, x, y, i, j, cnt, ix, iy;
    int x0, x1, y0, y1, y_cnt;
    int xs, xe, ys, ye;
    l32 al32Sum[IMAGE_WIDTH];
    float max_discrim;	

    GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;
	GTargetInfo  &TargetInfo = g_GblobalSettings.TargetInfo;
    TRect * ptRectRoi = &g_GblobalSettings.ROIInfo.tRectRoi;      //ROI外接矩形

    TRect atAllRect[128];
    float afScore[128];
    l32 l32RectCount = 0;
    l32 l32TargetNum;

    f32 f32DetectArrayNormW = VanishInfo.f32HumanNormW;
    f32 f32DetectArrayNormH = VanishInfo.f32HumanNormH;

    // 判断栅格底边的扫描起始位置
    f32 f32MinYs1 = VanishInfo.GetVanishY0ByH(f32DetectArrayNormH, 6);
    f32 f32MinYs2 = VanishInfo.GetVanishY0ByW(f32DetectArrayNormW, 10);
    ys = MAX(f32MinYs1, f32MinYs2);
    if(ys < ptRectRoi->l32Top)  ys = ptRectRoi->l32Top;
 
    // 判断栅格底边的扫描结束位置
    ye = MIN(IMAGE_HEIGHT, ptRectRoi->l32Top + ptRectRoi->l32Height);

    // 	UmpPrintf("ys = %d, ye = %d, g_fVanishY0 = %f, l32Top = %d, height = %d\n", ys, ye, g_fVanishY0,
    // 		ptTrackKLTHandle->ptRectRoi->l32Top, ptTrackKLTHandle->ptRectRoi->l32Height);

    //限制帧间差的绝对值, 防止过大的减除值影响后继判断
    for(int i = 0; i < IMAGE_HEIGHT * IMAGE_WIDTH; i++)
    {
        pu8Diff[i] = pu8Diff[i] > 50 ? 50 : pu8Diff[i];
    }

	//画框目标大小一致时，栅栏法失效
	if(ys >= ye)
	{
		return;
	}

#if OPENCV_DEBUG == 1
    Mat mImage;
    mImage = Mat::zeros(IMAGE_HEIGHT, IMAGE_WIDTH, CV_8UC1);
    memcpy(mImage.data, pu8Diff, IMAGE_WIDTH * IMAGE_HEIGHT);
    //cvtColor(mImage, mImage, CV_GRAY2RGB);
#endif

    //感兴趣区域的（帧间差的绝对值）积分图像
    getII(g_au32IIimage + (ys * (IMAGE_WIDTH + 1)), pu8Diff + (ys * IMAGE_WIDTH),
        IMAGE_WIDTH, IMAGE_WIDTH, (ye-ys));

    u32 * pu32IILine0;
    u32 * pu32IILine1;
    //收集栅栏运动物体存在证据
    y_cnt = 0;
    for(y= ys; y<ye; )
    {
        cnt = 0;
        y1 = y;
        w = VanishInfo.GetW(f32DetectArrayNormW, y1);
        h = VanishInfo.GetH(f32DetectArrayNormH, y1);
		
        y0 = y - h;

        if(y0 < ys)
        {
            y++;
            continue;
        }

        max_discrim = 0;
        xs = 0;
        xe = IMAGE_WIDTH - w;
        pu32IILine0 = g_au32IIimage + (IMAGE_WIDTH + 1)*(y0);
        pu32IILine1 = g_au32IIimage + (IMAGE_WIDTH + 1)*(y1);
        for(x=xs; x<xe; x+=w/2)
        {
            al32Sum[cnt] = 0;
            x0 = x;
            x1 = x + w;
            al32Sum[cnt] = pu32IILine0[x0] + pu32IILine1[x1] - pu32IILine0[x1] - pu32IILine1[x0];
            //printf("%d,", al32Sum[cnt]);
            cnt++;
        }
        //printf("\r\n");
        for(x=xs,i = 0; x<xe; x+=w/2,i++)
        {
            TRect tObjRc;
            tObjRc.l32Left = x;
            tObjRc.l32Top = y0;
            tObjRc.l32Width = w;
            tObjRc.l32Height = h;

            //排除当前已经跟踪的目标
            if(GetTRectOverlapArea(&tObjRc, pCurRect) > 0)
			{
				continue;
			}

			//排除目标为车的区域
			l32TargetNum = 0;
			for(j = 0; j < TargetInfo.l32TargetNum; j++)
			{
				if(GetTRectOverlapArea(&tObjRc, &TargetInfo.atRect[j]) > 0)
				{
					l32TargetNum++;
				}
			}

			if(l32TargetNum > 0)
			{
				continue;
			}

            int l32RoiNum;
            if(IsInROI(x + w/2, y0 + h/2, &l32RoiNum) == -1) 
				continue;

            float discrim = 0;
            if(i > 0 && i<cnt-1)
            {
                int da = al32Sum[i] - (al32Sum[i-1]);
                int db = al32Sum[i] - (al32Sum[i+1]);
                discrim = da < db?da:db;
            }
            discrim = discrim/(w*h);
#define DISCRIM_TH  3
#if OPENCV_DEBUG == 1
            Rect rc;
            rc.x = x;
            rc.y = y - h;
            rc.width = w;
            rc.height = h;
            if(discrim > DISCRIM_TH)
            {
                //printf("(%d,%d,%d)", al32Sum[i-1],al32Sum[i],al32Sum[i+1]);
                // 				 char strText[256];
                // 				 sprintf(strText, "%d", (int)discrim);
                // 				 putText(mImage, strText, Point(x, y), FONT_HERSHEY_SIMPLEX,0.4f,CV_RGB(255,0,0));
                rectangle(mImage, rc, CV_RGB(0,0,0));
            }
            else
            {
                // 				 char strText[256];
                // 				 sprintf(strText, "%d", (int)discrim);
                // 				 putText(mImage, strText, Point(x, y-5), FONT_HERSHEY_SIMPLEX,0.4f,CV_RGB(0,0,255));
                rectangle(mImage, rc, CV_RGB(0,0,255));
            }
#endif

            //对于本行栅栏中发现的物体，选择一个最有区分度的返回
            if(discrim > DISCRIM_TH)
            {
                int leftShoulderOffset,rightShoulderOffset;
                findShoulderBorder(pu8Diff, IMAGE_WIDTH, IMAGE_HEIGHT,
                    x, y-h, x+w, y,
                    leftShoulderOffset,rightShoulderOffset);

                float fOffsetRatio = (float)MAX(leftShoulderOffset , rightShoulderOffset)/(w);

                /*
                if(*pl32RectCnt < MAX_FIND_TARGET_NUM)
                {
                pRect[*pl32RectCnt] = tObjRc;
                afScore[*pl32RectCnt] = fOffsetRatio;
                (*pl32RectCnt) ++;
                }
                */

                if(fOffsetRatio < 1.0 && (max_discrim < discrim))				
                {
					//找到有效目标区域时，判断目标区域图像的平坦性，确定是否为真实目标
                    max_discrim = discrim;
                    if(l32RectCount < (sizeof(atAllRect)/sizeof(atAllRect[0])))
                    {
                        atAllRect[l32RectCount] = tObjRc;
                        afScore[l32RectCount] = fOffsetRatio;
						l32RectCount ++;
                    }
                }
            }
        }
        y += h;
    }

#if OPENCV_DEBUG == 1
    //imshow("GetFPSobj_vehicle", mImage);
    ShowY8InOpenCV("detectArray", mImage.data,IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH, 1);
    ShowRgbRec3("DetectMovingHumanByFENCE", pu8Diff, IMAGE_WIDTH, IMAGE_HEIGHT, 1, l32RectCount, atAllRect, afScore);
#endif

    if(l32RectCount > 0)
    {
        int id = rand() % l32RectCount;
        pNewRect[0] = atAllRect[id];
        printf("      fOffsetRatio = %f\n", afScore[id]);
    }
    /*
    // 按区域挑选
    ptTrackKLTHandle->l32NewRectCnt = 0;
    for(int i = 0; i < (*pl32RectCnt); i++)
    {
    if(pRect[i].l32Width > 0 && i > ptTrackKLTHandle->l32SelectRectYNo)
    {
    memcpy((void *)&ptTrackKLTHandle->tNewRect, &pRect[i], sizeof(TRect));
    ptTrackKLTHandle->l32NewRectCnt = 1;
    ptTrackKLTHandle->l32SelectRectYNo = i;
    break;
    }
    }

    if(ptTrackKLTHandle->l32NewRectCnt == 0)
    {
    for(int i = 0; i < (*pl32RectCnt); i++)
    {
    if(pRect[i].l32Width > 0)
    {
    memcpy((void *)&ptTrackKLTHandle->tNewRect, &pRect[i], sizeof(TRect));
    ptTrackKLTHandle->l32NewRectCnt = 1;
    ptTrackKLTHandle->l32SelectRectYNo = i;
    break;
    }
    }
    }
    */
    return;
}

static void RectExpand(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs)
{
    l32 al32Sum[IMAGE_WIDTH];
    int xs, xe;
    memcpy((void *)pNewRect, (void *)pOldRect, sizeof(TRect));

    GSVanishInfo &VanishInfo = g_GblobalSettings.VanishInfo;

    //float ratio = (pOldRect->l32Top +pOldRect->l32Height - VanishInfo.f32Y0)/(VanishInfo.f32TargetY - VanishInfo.f32Y0);
    float fvehicle_w = VanishInfo.GetW(VanishInfo.f32SaloonNormW, pOldRect->l32Top +pOldRect->l32Height); //VanishInfo.f32MaxW * ratio;
    float fvehicle_h = VanishInfo.GetH(VanishInfo.f32SaloonNormH, pOldRect->l32Top +pOldRect->l32Height);

    // 计算向左扩出的w/2每一列的最值投影
    int expand_max = fvehicle_w - pOldRect->l32Width;
    if(expand_max < 0) expand_max = pOldRect->l32Width/2;

    xs = pOldRect->l32Left - expand_max;
    xs = xs < 0 ? 0 : xs;
    xe = pOldRect->l32Left;
    for(int x=xe-1; x>=xs; x--)
    {
        u8 u8GrayMax = 0;
        for(int y = pOldRect->l32Top; y < pOldRect->l32Top + pOldRect->l32Height; y++)
        {
            if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
            {
                u8GrayMax = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
        if(u8GrayMax > 50)
        {
            pNewRect->l32Left--;
            pNewRect->l32Width++;
            expand_max --;
        }
        else
        {
            break;
        }
    }
    xs = pOldRect->l32Left + pOldRect->l32Width;
    xe = xs + expand_max;
    xe = xe > (IMAGE_WIDTH-1)?(IMAGE_WIDTH-1):xe;
    for(int x=xs; x<xe; x++)
    {
        u8 u8GrayMax = 0;
        for(int y = pOldRect->l32Top; y < pOldRect->l32Top + pOldRect->l32Height; y++)
        {
            if(pu8DiffAbs[y * IMAGE_WIDTH + x] > u8GrayMax)
            {
                u8GrayMax = pu8DiffAbs[y * IMAGE_WIDTH + x];
            }
        }
        if(u8GrayMax > 50)
        {
            pNewRect->l32Width++;
        }
        else
        {
            break;
        }
    }


}

void RectShrink(TRect *pOldRect, TRect *pNewRect, u8 *pu8DiffAbs)
{
    TRect atRect[2];
    int x, x1, x0, y, y0, y1, sum;
    u32 *pu32IILine0, *pu32IILine1;
    l32 al32Sum[IMAGE_WIDTH];

    y0 = pOldRect[0].l32Top;
    y1 = pOldRect[0].l32Top + pOldRect[0].l32Height;
    for(x=0; x<pOldRect[0].l32Width; x++)
    {
        x0 = x + pOldRect[0].l32Left;
        x1 = x0 + 1;
        sum = 0;
        u8 *pu8Temp = pu8DiffAbs + y0*IMAGE_WIDTH + x0;
        for(y=y0; y<y1; y++ )
        {
            sum += pu8Temp[0];
            pu8Temp += IMAGE_WIDTH;
        }
        al32Sum[x] = sum;
    }

    memcpy((void *)pNewRect, (void *)pOldRect, sizeof(TRect));
    for(x=0; x<pOldRect[0].l32Width; x++)
    {
        if(al32Sum[x] < 3*pOldRect->l32Height)
        {
            pNewRect->l32Left++;
            pNewRect->l32Width--;
        }
        else
        {
            break;
        }
    }
    for(x=pOldRect[0].l32Width - 1; x >= 0; x--)
    {
        if(al32Sum[x] < 3*pOldRect->l32Height)
        {
            pNewRect->l32Width--;
        }
        else
        {
            break;
        }
    }
}

void sobel_sum_hv(u8 *pu8Src, l32 l32Width, l32 l32Height, l32 l32Stride,
                  l32 &l32SobelPowerV, l32 &l32SobelPowerH)
{
    l32 l32OffsetU, l32OffsetD;
    l32 l32OffsetL, l32OffsetR;
    l32 l32I, l32J;
    l32 l32TempH, l32TempV;

    l32SobelPowerV = l32SobelPowerH = 0;

    for(l32I = 1; l32I < l32Height - 1; l32I++)
    {
        for(l32J = 1; l32J < l32Width - 1; l32J++)
        {
            l32OffsetU = (l32I - 1) * l32Stride + (l32J);
            l32OffsetD = (l32I + 1) * l32Stride + (l32J);
            l32TempH = pu8Src[l32OffsetU] + 2 * pu8Src[l32OffsetU + 1] + pu8Src[l32OffsetU + 2] -
                pu8Src[l32OffsetD] - 2 * pu8Src[l32OffsetD + 1] - pu8Src[l32OffsetD + 2];
            l32TempH = abs(l32TempH);

            l32OffsetL = (l32I) * l32Stride + (l32J - 1);
            l32OffsetR = (l32I) * l32Stride + (l32J + 1);
            l32TempV = pu8Src[l32OffsetL] - pu8Src[l32OffsetR] + 2 * (pu8Src[l32OffsetL + l32Stride]
            - pu8Src[l32OffsetR + l32Stride]) + pu8Src[l32OffsetL + 2 * l32Stride] - pu8Src[l32OffsetR + 2 * l32Stride];
            l32TempV = abs(l32TempV);


            l32SobelPowerH += l32TempH;
            l32SobelPowerV += l32TempV;
            /*
            if(l32TempH > l32TempV)
            {
            l32SobelPowerH += l32TempH + l32TempV;
            }
            else
            {
            l32SobelPowerV += l32TempH + l32TempV;
            }
            */
        }
    }
}






#endif
