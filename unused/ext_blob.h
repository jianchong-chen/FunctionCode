#ifndef ALTAIR_ALGM_H
#define ALTAIR_ALGM_H

/*#include <opencv2/opencv.hpp>*/

#define IMAGE_WIDTH     704
#define IMAGE_HEIGHT    576

#define MAX_NUM_BLOBS   128

#define	MAX_CC_IDS		100000	//2048 modified by  LDM

typedef struct _RunLength_
{
    int seg_id; //前景区间的pass1 id
    int seg_st; //前景区间的起始位置，指向前景区间的第一个前景位置
    int seg_ed; //前景区间的结束位置，指向前景区间之后的第一个背景位置
}TRunLength;

typedef struct AntBlob
{
   //POINT centroid;
   //POINT bbTopLeft, bbBottomRight;
    int x0;
    int x1;
    int y0;
    int y1;

    int mass;
    int id;           // connected component id
    int user;         // meaningless data storage for external use
    int isUse;
} AntBlob;

class CAntImage
{
public:
	CAntImage(int vw, int vh);
	~CAntImage();


    void ExtractBlobByMophCC(u8 * pu8FgMask);

	void setParameters(int targetSizeX,int targetSizeY, int xDilateSize, int yDilateSize); //set parameters of the image object

	//释放所有的内存
	void CAntImage::FreeAllMemory();
		
	int m_minTargetSizeX,m_minTargetSizeY;
	int m_x_dilate_size, m_y_dilate_size;
	unsigned char *m_pFgMask; // foreground mask
	unsigned char *m_pTempMask; // temp mask

public:
	//data members
	const int m_width; //image width
	const int m_height; //image height
	const int m_nChannels;

	int *m_pTempf1;  // 32 x 1    temp images

    AntBlob *m_blobs[MAX_NUM_BLOBS]; //blobs
	int m_nBlobs;
	int m_map[MAX_CC_IDS]; //map of connected components

	//private member functions
	BOOL erodeMask();
	BOOL dilateMask();
	BOOL findConnectedComponents();
	BOOL findBlobs();
    BOOL findBlobs2(TRunLength *ptLineRL, int nCCNumber);
};

#endif