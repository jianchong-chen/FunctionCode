


#define IMAGE_WIDTH 704
#define IMAGE_HEIGHT 576
///////////////////////////////////////////////////////////////////////////////////
//外部定义
//LCGMM参数
typedef struct tagLcgmmPara
{
    float f32AlphaPara;
    float f32RhoPara;
    float f32SigmaPara;
    float f32ThreshPara;
    float f32PixelMaxVar;
    float f32PixelMinVar;
}TLcgmmParam;

typedef struct
{
    float f32GmmPlusTs;
    float f32GmmPlusTg;
}TGmplusParam;

typedef struct
{
    TLcgmmParam tStLcg;
    TGmplusParam tStGmplus;
}TLcgmmPlusParam;

///////////////////////////////////////////////////////////////////////////////////
//内部定义
#define K_MODELS           3
#define bgT                0.79f  //背景权重阈值
#define InitialWeight      0.02f
#define InitialVariance    0.01f


typedef struct PixelModel
{
	float weight;
	float mean;
	float var;   
	float dist2;   
	float mah2;    

} PixelModel;


class CGMMPixel
{
public:
	CGMMPixel();
	~CGMMPixel();
	u8  fg;    
	void Process(float data, TLcgmmParam &Param);  
    void Process2(float data, TLcgmmParam &Param);  
public:
	PixelModel model[K_MODELS];
	PixelModel *pmMatch;
};


class CLcgmm
{
public:
    CLcgmm(int l32Width=IMAGE_WIDTH, int l32Height=IMAGE_HEIGHT);
    ~CLcgmm();

    void process(u8 *pu8Img, u8 *pu8MeanBg, u8 *pu8Mask);

    void setParameters(TLcgmmPlusParam *ptLcgmmPara);

private:
    void lcgmm(u8 *pu8Y);
    void gmmplus(u8 *pu8YUV, u8 *pu8MeanBg, u8 *pu8Mask);
    void UpdateBgImg(u8 *pu8YUV, u8 *pu8MeanBg, u8 *pu8Mask);

public:
    //data members
    const int m_l32Width; //image width
    const int m_l32Height; //image height

    float m_f32LearnRate; // Gmm+背景的学习率
    float m_f32LcLearnRate; // 缩小图像的Gmm+背景的学习率
    
    float m_f32Ts;
    float m_f32Tg;

private:
    TLcgmmPlusParam m_tLcgmmPara;
    s16 *m_ps16Iimg;
    u32 *m_pu32Iimg;
    CGMMPixel **m_ppcGmmpixel;

    l32 m_l32BlockWidth;
    l32 m_l32BlockHeight;
    l32 m_l32StepX;
    l32 m_l32StepY;
    l32 m_l32NumBlocksH;
    l32 m_l32NumBlocksW;

    s8 *m_ps8GmmBuf;
    s16 *m_ps16VData;
    s16 *m_ps16HData;
    float *m_pfMeanBg; // 平均背景
};
