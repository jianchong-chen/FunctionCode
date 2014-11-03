// framePCA.cpp : 定义控制台应用程序的入口点。
#include "stdafx.h"

#include "LHVideoWriter.h"
#include "LHIPCA4BGModel.h"
#include <stdio.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <wchar.h>
using namespace std;

//#define GAUSSIAN_SMOOTH 13

int _tmain(int argc, char* argv[])
{
	__int64 t1, t2;
	float tF = cvGetTickFrequency();
	float avg_fps = 0;
	float tSum1 = 18;
	float tSum2 = 0;
	float tSum3 = 0;

	char buf_paras[10];
	char buf_filename[250];
	char save_filename[150];
	char videoDir[400];
	char resultDir[400];
	char videoName[100];
	char *captureSrc = NULL;
	char *resultVideo = NULL;	
	char ch;

	int frameCount = 0;
	CvCapture *capture = 0;
	//IplImage *frame_image = 0;
	IplImage *frame_image_gray = NULL;
	IplImage *mask_image = NULL;
	IplImage* reconError_image = NULL;
	IplImage *reconBGImage = NULL;
	IplImage *err_image = NULL;
	IplImage *large_image = NULL;
	IplImage *fg_objects = NULL;
	IplImage *ROI_Image = NULL;

	CLHVideoWriter videoWriter;
	int codeType = 0;

AA:	cout<<"Begin?(y/n): ";
	cin>>ch;
	if(ch != 'y')
	{
		goto AA;
	}

	int LOOP = 5;

	for ( int i = 0; i< LOOP; i++)
	{
		strcpy_s(save_filename, "IPCAResults_443_Sub2_Init50_r200_Adap10_min16_Alpha1_PC5_");
		
		ifstream fVideoList;
		fVideoList.open( "E:\\SABS\\List.txt", ios::out );

		char* brick_scale = "4X4X3"; 
		for ( int nv = 0; nv < 5; nv++)
		{
			cout<<"Video:"<<nv+1<<"******************************\n";

			strcpy_s(buf_filename, save_filename);
			strcpy_s(videoDir, "E:\\SABS\\");
			strcpy_s(resultDir, "E:\\SABS\\Results\\");

			fVideoList.getline( videoName, 99, '\n');
			strcat_s( videoDir, videoName );
			strcat_s( buf_filename, videoName );
			strcat_s( resultDir, buf_filename );
			captureSrc = videoDir;
			resultVideo = resultDir;

			// Capture video frames
			capture = 0;
			if(strlen(captureSrc) == 1 && isdigit(captureSrc[0]))
			{
				capture = cvCaptureFromCAM( captureSrc[0]-'0');
			}
			else 
			{
				capture = cvCaptureFromAVI( captureSrc ); 
			}
			if(!capture)
			{
				fprintf(stderr, "Could not initialize capturing...\n");
				return false;
			}

			IplImage *frame_image = 0;		
			frame_image = cvQueryFrame(capture);    // not be Released!?
			if(!frame_image)
			{
				printf("Video Read Error!\n");
				continue;
			}

			frame_image_gray = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 1);
			cvCvtColor(frame_image, frame_image_gray, CV_RGB2GRAY);

			mask_image = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 1);
			reconError_image = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 1);
			reconBGImage = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 1);
			err_image = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 1);
			large_image = cvCreateImage(cvSize(frame_image->width*3, frame_image->height*2), 8, 3);
			fg_objects = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 3);		

			//codeType = CV_FOURCC('X','V','I','D');
			codeType = CV_FOURCC('D','I','V','X');
			int codeType_2 = cvGetCaptureProperty( capture, CV_CAP_PROP_FOURCC );
			videoWriter.Init( resultVideo, codeType, 25, frame_image->height*2, frame_image->width*3, 1 );

			ROI_Image = cvLoadImage( "E:\\roi_image.bmp", 0 );
			if(ROI_Image == 0)
			{
				ROI_Image = cvCreateImage(cvSize(frame_image->width, frame_image->height), 8, 1);
				for(int i=0; i<ROI_Image->height; i++)
				{
					for(int j=0; j<ROI_Image->width; j++)
					{
						CV_IMAGE_ELEM(ROI_Image, uchar, i, j) = 255;
					}
				}
			}

			int nInitialSample = 50;     // Number of frames used for initializing
			int updateSpeed = -1;        //正常更新速度每帧7行,为-1时意味需要新一次更一帧所有行
			CLHIPCA4BGModel onlinePCAModel;
	
			onlinePCAModel.Initial(frame_image_gray, ROI_Image, nInitialSample);

			CvFont font;
			cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 0.7, 0.7, 0, 1.5);
			CvFont pFont;
			cvInitFont(&pFont,CV_FONT_HERSHEY_PLAIN,2.0,2.0,0,2);
			char buf_int[15];

			frameCount = 0;
			while(frame_image = cvQueryFrame(capture))
			{ 
				
				frameCount++;
				//cout<<"frameCount: "<<frameCount<<endl;		
				cvFlip(frame_image);
				cvCvtColor(frame_image, frame_image_gray, CV_RGB2GRAY);
#ifdef GAUSSIAN_SMOOTH
				cvSmooth(frame_image_gray, frame_image_gray, CV_GAUSSIAN, GAUSSIAN_SMOOTH, GAUSSIAN_SMOOTH);
#endif

				t1 = cvGetTickCount();
				cvZero(mask_image);
				cvSetZero(reconError_image);
				if(onlinePCAModel.foregroundDetect(frame_image_gray, mask_image, reconError_image, updateSpeed))
				{
					t2 = cvGetTickCount();
					//cout<<"Loop(): "<<(t2-t1)/(tF*1000)<<"ms"<<endl;
					if(frameCount == nInitialSample+10)
					{
						tSum1 = (t2-t1)/(tF*1000.0);
					}
					else if(frameCount > nInitialSample+10)
					{
						//tSum1 = (1-1/(frameCount-200.0))*tSum1 + ((t2-t1)/(tF*1000.0))/(frameCount-200.0);
						tSum1 = (1-0.05)*tSum1 + 0.05*(t2-t1)/(tF*1000.0);
						avg_fps = 1000.0/tSum1;
					}	

					cvSetZero(large_image);
					cvSetZero(fg_objects);
					cvSetZero(reconBGImage);
					cvSetZero(err_image);				

					if (frameCount > nInitialSample+10)
					{
						onlinePCAModel.reconBGFrame(reconBGImage);
					}

					cvSetImageROI( large_image, cvRect(0,0,frame_image->width,frame_image->height) );
					cvCopy(frame_image,large_image);
					cvResetImageROI(large_image);
					cvSetImageROI( large_image, cvRect(frame_image->width,0,frame_image->width,frame_image->height) );
					cvCvtColor( mask_image, large_image, CV_GRAY2RGB);
					cvResetImageROI(large_image);
					cvSetImageROI( large_image, cvRect(0,frame_image->height,frame_image->width,frame_image->height) );
					cvCopy(frame_image,fg_objects,mask_image);
					cvResize(fg_objects,large_image);
					cvResetImageROI(large_image);
					cvSetImageROI( large_image, cvRect(frame_image->width,frame_image->height,frame_image->width,frame_image->height) );
					cvCvtColor( reconBGImage, large_image, CV_GRAY2RGB);
					cvResetImageROI(large_image);
					cvSetImageROI( large_image, cvRect(frame_image->width*2, 0, frame_image->width, frame_image->height) );
					cvCvtColor( reconError_image, large_image, CV_GRAY2RGB);
					cvResetImageROI(large_image);

					//cvNamedWindow("Error",0);
					//cvShowImage("Error",reconError_image);

					_itoa_s(frameCount,buf_int,10);
					cvPutText(large_image, buf_int, cvPoint(10,30), &pFont, CV_RGB(255,0,0));

					//_itoa(matchThreshold,buf_int,10);
					cvPutText(large_image, /*buf_int*/brick_scale, cvPoint(100,30), &pFont, CV_RGB(0,255,0));

					//_itoa(avg_fps,buf_int,10);
					sprintf_s(buf_int, "%.2f", avg_fps); 
					cvPutText(large_image, buf_int, cvPoint(210,30), &pFont, CV_RGB(255,0,0));

					cvPutText(large_image, "fps", cvPoint(285,30), &pFont, cvScalar(0,255,255));

					videoWriter.WriteFrame(large_image->imageData);

					cvNamedWindow("Result",0);
					cvShowImage("Result",large_image);
					cvWaitKey(1);
				}
			}

			videoWriter.Release();		
			onlinePCAModel.Release();

			cvReleaseImage(&frame_image_gray);
			cvReleaseCapture(&capture);
			cvReleaseImage(&ROI_Image);
			cvReleaseImage(&mask_image);
			cvReleaseImage(&reconError_image);
			cvReleaseImage(&err_image);
			cvReleaseImage(&large_image);
			cvReleaseImage(&fg_objects);
		}	

		fVideoList.close();
	}
	
	return 0;
}
