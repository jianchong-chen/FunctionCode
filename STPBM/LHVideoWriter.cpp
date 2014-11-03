/** LHVideoWriter.cpp
*
* Copyright (C) 2006~2007 Lotushill Institute
*
* Created: Jan. 9, 2007 by Kylin
* Modified:
*
* This file contains the implementation of VideoWriter Class of the Intelligent 
* Survaillance System.
*
* Implementation the VideoWriter method
* 
* @see 
* @author Kylin 
* @email lqi.lhi@gmail.com
*/
#include "stdafx.h"
#include "LHVideoWriter.h"


/////////////////////////////// PUBLIC ///////////////////////////////////////
//============================= LIFECYCLE ====================================
/** Default Constructor
*/

CLHVideoWriter::CLHVideoWriter()
{
}

/*virtual*/ CLHVideoWriter::~CLHVideoWriter()
{
}

//============================= OPERATORS ====================================
/* No Operators*/

//============================= OPERATIONS ===================================

/** Initialize
*
* PRECONDITION
* REQUIRE(isColor = 0 / 1)
*
* WARNING
* EXAMPLES
*
* @param filename       the name of the output video file
* @param fps            frame per second
* @param height         image height
* @param width          image width
* @param isColor        whether the image is color
*
* @return void
*/
void CLHVideoWriter::Init(const char* filename, int codeType, float fps, int height, int width, int isColor)
{
    mHeight = height;
    mWidth = width;
    if (isColor != 0)
    {
        mpFrameImage = cvCreateImageHeader(cvSize(width, height), IPL_DEPTH_8U, 3);
    }
    else
    {
        mpFrameImage = cvCreateImageHeader(cvSize(width, height), IPL_DEPTH_8U, 1);
    }
    //mpVideoWriter = cvCreateVideoWriter(filename, CV_FOURCC('D','I','V','X'), fps, cvSize(width, height), isColor);
	//mpVideoWriter = cvCreateVideoWriter(filename, CV_FOURCC('X','V','I','D'), fps, cvSize(width, height), isColor);
	mpVideoWriter = cvCreateVideoWriter(filename, codeType, fps, cvSize(width, height), isColor);
}

/** Writer a frame image to video file
*
* PRECONDITION
* REQUIRE
*
* WARNING
* EXAMPLES
*
* @param imageData          a pointer holds the pixel data
*
* @return void
*/
void CLHVideoWriter::WriteFrame( char* imageData)
{
    mpFrameImage->imageData = (char*)imageData;
    cvWriteFrame(mpVideoWriter, mpFrameImage);
}

/** Release resources
*
* PRECONDITION
* REQUIRE
*
* WARNING
* EXAMPLES
*
* @param void
*
* @return void
*/
void CLHVideoWriter::Release()
{
    cvReleaseImageHeader(&mpFrameImage);
    cvReleaseVideoWriter(&mpVideoWriter);
}
