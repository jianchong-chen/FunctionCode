/** LHVideoWriter.h
*
* Copyright (C) 2006~2007 Lotushill Institute
*
* Created: Jan, 2006 by Kylin
* Modified:
*
* This class is for writing a frame image to a video file.
*
* #include "LHVideoWriter.h" <BR>
* Corresponding implementation file is "LHVideoWriter.cpp"
*
* @see 
* @author Kylin 
* @email lqi.lhi@gmail.com
*/

#ifndef LHVideoWriter_H
#define LHVideoWriter_H
#include "stdafx.h"
#include "highgui.h"

#ifdef CImage
#undef CImage
#endif

class CLHVideoWriter
{
public:
    /// Default constructor
    CLHVideoWriter(void);

    /// Destructor
    virtual ~CLHVideoWriter(void);

    /// Initialize
    void Init(const char* filename, int codeType, float fps, int height, int width, int isColor);

    /// Write a frame image to the video file
    void WriteFrame( char* imageData);

    /// Release resources
    void Release();

private:

    /// Pointer of video writer in OpenCV lib
    CvVideoWriter* mpVideoWriter;

    /// Image width
    int mWidth;

    /// Image height
    int mHeight;
    
    /// IplImage pointer of a frame image
    IplImage* mpFrameImage;
};

#endif
