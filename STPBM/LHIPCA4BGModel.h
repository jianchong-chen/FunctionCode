/** CLHIPCA4BGModel.h
*
* Copyright (C) 2007~2008 Lotushill Institute
*
* Created: Mar, 2008 by Youdong Zhao
* Modified:
*
* *******
* #include "LHIPCA4BGModel.h" <BR>
* Corresponding implementation file is "LHIPCA4BGModel.cpp"
*
*/

#pragma once
//#include "LHParas.h"
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "LHCCIPCA.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

class CLHIPCA4BGModel
{
public:

	// initialize the IPCABGModel using these parameters
	bool Initial( IplImage *current_image, IplImage* ROI_Image, int nInitialObs);	

	// get the foreground mask image corresponding to the current frame
	bool getFGMask( const IplImage* current_image, IplImage* mask_image, IplImage* reconError_image, const int updateSpeed );

	// using the current frame to update the IPCABGModel and obtain the corresponding foreground mask image
	bool foregroundDetect( const IplImage* current_image, IplImage* mask_image, IplImage* reconError_image, const int updateSpeed );

	bool reconBGFrame(IplImage* reconBgImage);

	bool Release();

private:
	// The frame brick as a LoopLinkTable
	struct Brick{
				int length;
				int ID;
				IplImage* Data;
				struct Brick* Next;
	}*pHeader;

	//static const int BRICK_STEP = 4;            // SubBrick's large size that used to compute
	//static const int PATCH_SIZE = 1;            // 提取Brick时的步宽，一定为1，别设大的，不然可能有问题，m_half_width_brick的计算与他有关。为1时可以保证没问题
	//static const int BRICK_PLY = 3;             // SubBrick's ply
	//static const int TARGET_DIM = 5;            // the target dimensions of onlinePCA subspace
	//static const int SAMPLE_SIZE = 2;           // Down-Sample rate

	int BRICK_STEP;           // SubBrick's large size that used to compute
	int PATCH_SIZE;           // 提取Brick时的步宽，一定为1，别设大的，不然可能有问题，m_half_width_brick的计算与他有关。为1时可以保证没问题
	int BRICK_PLY;            // SubBrick's ply
	int TARGET_DIM;           // the target dimensions of onlinePCA subspace
	int SAMPLE_SIZE;          // Down-Sample rate

	int m_nInitialObs;         // the number of samples that used to initialize the onlinePCA model
	float m_forgetFactor;	   // the forgotten factor that used by the onlinePCA_update algorithm
	float m_Threshold_Divid;      // sub_brick match threshold;
	int m_half_width_brick;

	__int64 frameCount;
	int loop_flag;
	bool matched_flag;

	int sizes_brick[3];       // Video frame brick's Size
	int subBrickDims;	     // Used dims of sub brick
	int orig_subBrickDims;   // Original dims of sub brick
	CvMat *sub_brick;       // Used sub brick     //用来解决brick的长度(dims)为非4的倍数问题，因为后面计算用了SSE加速计算，要求dims为4的倍数，即舍去几个像素，如27-3=24
	CvMat *orig_sub_brick;  // Original sub brick
	CvMat *recon_brick;

	float m_max_recon_error;
	CvMat *reconError_mat;
	CvMat *reconError_mat_255;
	IplImage* mask_image_temp;
	IplImage* gray_image;
	//IplImage* gray_image_small;
	IplImage* current_image_small;
	IplImage* recon_image_small;

	int MM;
	int NN;

	vector<CLHCCIPCA> pca_model;	
	vector<int> fgCount;

	CvMat* ROI_Mat;
	// obtain the sub_brick from the whole Brick w.r.t.the indexs ii and jj 
	bool ObtainSubBrick(const int ii, const int jj);

	// convert frames to a Brick
	bool frame2Brick(const IplImage* current_image);
};