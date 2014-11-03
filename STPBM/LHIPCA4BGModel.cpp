// framePCA.cpp : 定义控制台应用程序的入口点。
#include "stdafx.h"

#include "LHParas.h"
#include "LHIPCA4BGModel.h"
#include <fstream>
using namespace std;

bool CLHIPCA4BGModel::ObtainSubBrick( const int ii, const int jj )
{
	int ii_start = ii*PATCH_SIZE-(BRICK_STEP-PATCH_SIZE)/2;   // 提取Brick数据的起始点，由中心(ii，jj)向前减去一半宽度
	int jj_start = jj*PATCH_SIZE-(BRICK_STEP-PATCH_SIZE)/2;
	int row_step = ii_start+BRICK_STEP;
	int col_step = jj_start+BRICK_STEP;
	int cols = 0;
	struct Brick* pTemp = pHeader;
	int numChannels = pTemp->Data->nChannels;

	do{
		for(int i=ii_start; i<row_step; i++)
		{
			for(int j=jj_start; j<col_step; j++)
			{
				for(int k = 0; k<numChannels; k++)
				{
					float temp = ((uchar*)(pTemp->Data->imageData + i*pTemp->Data->widthStep))[j*numChannels + k];
					CV_MAT_ELEM(*orig_sub_brick, float, 0, cols) = temp;
					cols++;
				}
			}
		}
		pTemp = pTemp->Next;
	}while(pTemp!=pHeader);

	for (int i=0; i<subBrickDims; i++)
	{
		CV_MAT_ELEM(*sub_brick, float, 0, i) = CV_MAT_ELEM(*orig_sub_brick, float, 0, i);
	}
	return true;
}

bool CLHIPCA4BGModel::frame2Brick(const IplImage* current_image)
{
	frameCount++;
		
	// obtain current brick 
	cvResize(current_image, current_image_small, CV_INTER_AREA);

	if(frameCount < sizes_brick[2])
	{
		cvScale(current_image_small, pHeader->Data, 1);
		pHeader = pHeader->Next;
		return false;
	}
	else
	{
		cvScale(current_image_small, pHeader->Data, 1);
		pHeader = pHeader->Next;
		return true;
	}
}

bool CLHIPCA4BGModel::Initial( IplImage* current_image, IplImage* ROI_Image, int nInitialObs)
{
	BRICK_STEP = LH_BRICK_STEP;           // SubBrick's size that used to compute
	PATCH_SIZE = LH_PATCH_SIZE;           // SubBrick's step that used to sample
	BRICK_PLY = LH_BRICK_PLY;            // SubBrick's ply
	TARGET_DIM = LH_TARGET_DIM;          // the target dimensions of onlinePCA subspace
	SAMPLE_SIZE = LH_SAMPLE_SIZE;

	m_nInitialObs = nInitialObs;       // the number of samples that used to initialize the onlinePCA model	
	m_forgetFactor = 0;                // the forgotten factor that used by the onlinePCA_update algorithm
	
	m_Threshold_Divid = 0;
	
	m_max_recon_error = 0;
	m_half_width_brick = (BRICK_STEP - PATCH_SIZE)/2 + 1;

	frameCount = 0;
	loop_flag = 2;
	matched_flag = false;

	sizes_brick[0] = current_image->height/SAMPLE_SIZE;      // Video frame brick's Height
	sizes_brick[1] = current_image->width/SAMPLE_SIZE;       // Video frame brick's Width
	sizes_brick[2] = BRICK_PLY;				                 // Video frame brick's Ply

	orig_subBrickDims = BRICK_STEP*BRICK_STEP*BRICK_PLY*current_image->nChannels;
	subBrickDims = orig_subBrickDims - orig_subBrickDims % 4;
	orig_sub_brick = cvCreateMat(1, orig_subBrickDims, CV_32FC1);
	sub_brick = cvCreateMat(1, subBrickDims, CV_32FC1);
	recon_brick = cvCreateMat(sub_brick->height, sub_brick->width, CV_32FC1);

	MM = sizes_brick[0]/PATCH_SIZE;
	NN = sizes_brick[1]/PATCH_SIZE;

	for( int ii = 0; ii < MM; ii++ )
	{
		for( int jj = 0; jj < NN; jj++ )
		{  
			CLHCCIPCA tempModel;
			pca_model.push_back(tempModel);
		}
	}

	for( int ii = 0; ii < MM; ii++ )
	{
		for( int jj = 0; jj < NN; jj++ )
		{  
			int tempInt = 0;
			fgCount.push_back(tempInt);
		}
	}

	ROI_Mat = cvCreateMat(MM, NN, CV_8UC1);
	cvResize(ROI_Image, ROI_Mat, 1);
	for(int i=0; i<ROI_Mat->rows; i++)
	{
		for(int j=0; j<ROI_Mat->cols; j++)
		{
			if( CV_MAT_ELEM(*ROI_Mat, uchar, i, j) > 0 )
			{
				CV_MAT_ELEM(*ROI_Mat, uchar, i, j) = 255;
			}
		}
	}
	mask_image_temp = cvCreateImage(cvSize(NN, MM), current_image->depth, 1);
	reconError_mat = cvCreateMat( MM, NN, CV_32FC1);
	reconError_mat_255 = cvCreateMat( MM, NN, CV_8UC1);

	gray_image = cvCreateImage(cvSize(current_image->width, current_image->height), 8, 1);
	//gray_image_small = cvCreateImage(cvSize(sizes_brick[1], sizes_brick[0]), 8, 1);

	current_image_small = cvCreateImage(cvSize(sizes_brick[1], sizes_brick[0]), current_image->depth, current_image->nChannels);
	recon_image_small = cvCreateImage(cvSize(sizes_brick[1], sizes_brick[0]), current_image->depth, current_image->nChannels);

	// Initialize the brick as a LoopLinkTable
	struct Brick *pNode1, *pNode2; 
	pNode1 = pNode2 = new struct Brick;
	pNode1->ID = 1;
	pNode1->length = sizes_brick[2];
	pNode1->Data = cvCreateImage(cvSize(sizes_brick[1], sizes_brick[0]), current_image->depth, current_image->nChannels);
	for( int i=1; i<=sizes_brick[2]; i++)
	{
		if( i == 1 )
		{
			pHeader = pNode1;
		}
		else
		{
			pNode2->Next = pNode1;
		}

		pNode2 = pNode1;

		pNode1 = new struct Brick;
		pNode1->ID = i+1;
		pNode1->length = sizes_brick[2];
		pNode1->Data = cvCreateImage(cvSize(sizes_brick[1], sizes_brick[0]), current_image->depth, current_image->nChannels);
	}
	pNode2->Next = pHeader;

	cvReleaseImage(&pNode1->Data);
	delete pNode1;
	return true;
}

bool CLHIPCA4BGModel::Release()
{
	cvReleaseMat(&sub_brick);
	cvReleaseMat(&orig_sub_brick);
	cvReleaseMat(&recon_brick);
	cvReleaseImage(&current_image_small);
	cvReleaseImage(&recon_image_small);	
	cvReleaseImage(&gray_image);
	cvReleaseImage(&mask_image_temp);
	cvReleaseMat(&reconError_mat);
	cvReleaseMat(&reconError_mat_255);

	struct Brick* pTemp = pHeader;
	do{
		cvReleaseImage(&pTemp->Data);
		pTemp = pTemp->Next;
	}while(pTemp!=pHeader);

	delete pHeader;

	for(int i=m_half_width_brick;i<MM-m_half_width_brick;i++)
	{
		for(int j=m_half_width_brick;j<NN-m_half_width_brick;j++)
		{
			if(CV_MAT_ELEM(*ROI_Mat, uchar, i, j) == 255)
			{
				pca_model[j+i*NN].Release();
			}
		}
	}

	cvReleaseMat(&ROI_Mat);

	return true;
}


bool CLHIPCA4BGModel::getFGMask( const IplImage* current_image, IplImage* mask_image, IplImage* reconError_image, const int updateSpeed )
{	
	float errRecon = 0;
	float max_error = 300;
	float min_error = 1000000;
	float interval_error = 0;
	
	float updateErr=0;
	float r = LH_T_MUST_UPDATE;  //200
	int temp_updateSpeed = 0;

	frame2Brick(current_image);

	// Determine the line numbers that should be update by this frame
	if (updateSpeed == -1)
	{ // Update the whole frame
		temp_updateSpeed = MM;
	}
	else
	{ // Update partly
		temp_updateSpeed = updateSpeed;
	}
	if( loop_flag >= MM-m_half_width_brick )
	{
		loop_flag = m_half_width_brick;
	}
	int update_begin_line = loop_flag;
	loop_flag += temp_updateSpeed;

	cvSetZero(mask_image_temp);
	cvSetZero(reconError_mat);
	cvSetZero(reconError_mat_255);

	if(frameCount < sizes_brick[2])
	{
		cvResize(mask_image_temp,gray_image);
		if(mask_image->nChannels == 1)
		{
			cvCopy(gray_image,mask_image);
		}
		else if(mask_image->nChannels == 3)
		{
			cvCvtColor(gray_image, mask_image, CV_GRAY2RGB);
		}
		return false;
	}

	if(frameCount <= m_nInitialObs + sizes_brick[2])
	{
		cvResize(mask_image_temp,gray_image);
		if(mask_image->nChannels == 1)
		{
			cvCopy(gray_image,mask_image);
		}
		else if(mask_image->nChannels == 3)
		{
			cvCvtColor(gray_image, mask_image, CV_GRAY2RGB);
		}

		// Initializing background model
		for( int ii = m_half_width_brick; ii < MM-m_half_width_brick; ii++ )
		{
			for( int jj = m_half_width_brick; jj < NN-m_half_width_brick; jj++ )
			{  
				if(CV_MAT_ELEM(*ROI_Mat, uchar, ii, jj) == 255)
				{
					ObtainSubBrick( ii, jj ); // ii = 1:70, jj = 1:94
					//t2 = cvGetTickCount();
					//cout<<"ObtainSubBrick(): "<<(t2-t1)/(tF*1000)<<"ms"<<endl;

					if(frameCount-BRICK_PLY == 0)  //Initialization: applying the memory for the online-pca model
					{
						pca_model[ jj+ii*NN ].Initial(sub_brick, TARGET_DIM, m_nInitialObs, m_forgetFactor);
						pca_model[ jj+ii*NN ].UpdatePCA( sub_brick, updateErr );
					}
					else  //Initialization: computing the offline-pca model for initializing the online-pca model
					{ 
						//t1 = cvGetTickCount();
						pca_model[ jj+ii*NN ].UpdatePCA( sub_brick, updateErr );
						//t2 = cvGetTickCount();
						//cout<<"pca_model(): "<<(t2-t1)/(tF*1000)<<"ms"<<endl;
					}
				}
			}
		}
		return false;
	}


	//ofstream fPCA;

	//char buf_filename[100] = "E:\\Temp\\";
	//char buf_int[10];
	//_itoa(frameCount,buf_int,10);
	//strcat(buf_int,".txt");
	//char *filename = strcat(buf_filename,buf_int);

	//fPCA.open(filename, ios::trunc);
	////fPCA.open("I:\\match_test.txt", ios::trunc);

	for( int ii = m_half_width_brick; ii < MM-m_half_width_brick; ii++ )
	{
		for( int jj = m_half_width_brick; jj < NN-m_half_width_brick; jj++ )
		{ 
			if(CV_MAT_ELEM(*ROI_Mat, uchar, ii, jj) == 255)
			{			
				ObtainSubBrick( ii, jj );

				//int temp2 = fgCount[ jj+ii*NN ];

				//Matching
				int result_match = pca_model[ jj+ii*NN ].MatchSubBrick( sub_brick, errRecon, m_Threshold_Divid);
				if( !result_match )
				{ // Matching failure, foreground					
					fgCount[ jj+ii*NN ]++;
					CV_IMAGE_ELEM(mask_image_temp,uchar,ii,jj) = 255;
					if( fgCount[ jj+ii*NN ] > r)  
					{ // If a position is FG for "r" frames, use it to update the model
						pca_model[ jj+ii*NN ].UpdatePCA( sub_brick, updateErr );
					}
				}
				else if( result_match == 1)
				{// Matching successfully, background, and if reconError is not very small (e.g.,< 5),then update by it
					if(fgCount[ jj+ii*NN ] > 0)
					{ // 该位置前景计数减一，直到为0，为0后不再减，即这个计数的某种作用是说明当前帧是“连续”多少帧前景
						fgCount[ jj+ii*NN ]--;
					}

					if( update_begin_line <= ii && ii <= loop_flag - 1)
					{ // Only update part of all the models, indicated by the "updateSpeed"
						pca_model[ jj+ii*NN ].UpdatePCA( sub_brick, updateErr );
					}										
				}
				else if( result_match == 2)
				{
					if(fgCount[ jj+ii*NN ] > 0)
					{ // 该位置前景计数减一，直到为0，为0后不再减，即这个计数的某种作用是说明当前帧是“连续”多少帧前景
						fgCount[ jj+ii*NN ]--;
					}
					//if( update_begin_line <= ii && ii <= loop_flag - 1)
					//{ // Only update part of all the models, indicated by the "updateSpeed"
					//	pca_model[ jj+ii*NN ].UpdatePCA( sub_brick, updateErr );
					//}
				}

				CV_MAT_ELEM(*reconError_mat, float, ii, jj) = errRecon;
				if ( max_error < errRecon)
				{
					max_error = errRecon;
				}
				if (min_error > errRecon)
				{
					min_error = errRecon;
				}
				//if (frameCount > 3000 && frameCount < 3500)
				//{
					//fPCA<<errRecon<<" ";
				//}				
			}			
		}
		//if (frameCount > 3000 && frameCount < 3500)
		//{
			//fPCA<<endl;
		//}
	}
	//fPCA.close();


	//if ( m_max_recon_error < max_error )
	//{
	//	m_max_recon_error = max_error;
	//}

	//char sNum[6];	
	//_itoa(frameCount, sNum, 10);
	//string pathname = "E:\\Temp\\";
	//string fileName(sNum);
	//pathname = pathname + fileName + ".mat";
	//cvSave(pathname.data(), reconError_mat);

	interval_error = max_error - min_error;	
	cvSubS( reconError_mat, cvScalarAll(min_error), reconError_mat );
	cvAddWeighted( reconError_mat, 255.0/interval_error, reconError_mat, 0.0, 0, reconError_mat );
	for( int ii = m_half_width_brick; ii < MM-m_half_width_brick; ii++ )
	{
		for( int jj = m_half_width_brick; jj < NN-m_half_width_brick; jj++ )
		{ 
			int temp = (uchar)CV_MAT_ELEM(*reconError_mat, float, ii, jj);	
			CV_MAT_ELEM(*reconError_mat_255, uchar, ii, jj) = temp;
		}
	}
	cvResize( reconError_mat_255, reconError_image, 1 );
	//cvNamedWindow("Error",0);
	//cvShowImage("Error",reconError_image);
	
	cvResize(mask_image_temp,gray_image,CV_INTER_NN);
	if(mask_image->nChannels == 1)
	{
		cvCopy(gray_image,mask_image);
	}
	else if(mask_image->nChannels == 3)
	{
		cvCvtColor(gray_image, mask_image, CV_GRAY2RGB);
	}		

	return true;
}

bool CLHIPCA4BGModel::foregroundDetect( const IplImage* current_image, IplImage* mask_image, IplImage* reconError_image, const int updateSpeed )
{		
	__int64 t1, t2;
	float tF = cvGetTickFrequency(); 

	t1 = cvGetTickCount();
	getFGMask( current_image, mask_image, reconError_image, updateSpeed );
	t2 = cvGetTickCount();

	return true;
}


bool CLHIPCA4BGModel::reconBGFrame( IplImage* reconBgImage )
{
	__int64 t1, t2;
	float tF = cvGetTickFrequency(); 

	// 443 --> 37 current pixel
	
	t1 = cvGetTickCount();
	//fPCA.open(filename, ios::trunc);
	//fPCA.open("I:\\update_mean.txt", ios::trunc);
	//fPCA2.open(filename2, ios::trunc);

	cvSetZero(recon_image_small);
	for( int ii = m_half_width_brick; ii < MM-m_half_width_brick; ii++ )
	{
		for( int jj = m_half_width_brick; jj < NN-m_half_width_brick; jj++ )
		{ 
			if(CV_MAT_ELEM(*ROI_Mat, uchar, ii, jj) == 255)
			{
				ObtainSubBrick( ii, jj );
				//reconstruction
				cvSetZero(recon_brick);
				pca_model[ jj+ii*NN ].ReconBrick(sub_brick, recon_brick);
				CV_IMAGE_ELEM(recon_image_small, uchar, ii, jj) = (uchar)CV_MAT_ELEM(*recon_brick, float, 0, 5);
				//fPCA<<updateErr<<" ";
				//fPCA2<<updateErrRate<<" ";
			}
		}
		//fPCA<<endl;
		//fPCA2<<endl;
	}
	//fPCA.close();
	//fPCA2.close();
	cvResize(recon_image_small, reconBgImage, CV_INTER_NN);
	t2 = cvGetTickCount();
	//cout<<"ReconBGFrame(): "<<(t2-t1)/(tF*1000)<<"ms"<<endl;
	return true;
}
