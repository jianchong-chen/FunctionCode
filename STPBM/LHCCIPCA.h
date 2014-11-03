/** LHCCIPCA.h
*
* Copyright (C) 2007~2008 Lotushill Institute
*
* Created: Jan, 2007 by Wenze Hu
* Modified: Mar, 2008 by Youdong Zhao
*
* This class is for carLight detect, and is suitable for general purpose ONLINEPCA.
* It is a modified implementaion of TPAMI paper "Candid Covariacne free Incremental PCA"
* #include "LHCCIPCA.h" <BR>
* Corresponding implementation file is "LHCCIPCA.cpp"
*
*/

#pragma once
#include <cxcore.h>
#include <cv.h>

class CLHCCIPCA
{

public:
	/* allocate memory 
	*  pFirstObs, first observation, should be a column vecotr,
	*  targetDim, to what lower dimention,
	*  batchFrame, using how many obs for initial the pca
	*/

	void Initial(CvMat* pFirstObs, int targetDim, int batchFrame, float forgetFattor ); 

	void Initial2(void);

	// update the background model (subspace)
	bool UpdatePCA( const CvMat* pCurrentObs, float &updateErr );

	// match the current observations to the pca model
	int MatchSubBrick( const CvMat* pCurrentObs, float &errRecon, float m_Threshold_Divid );

	// reconstruct the observation brick using the BG subspace
	bool ReconBrick( const CvMat* pCurrentObs, CvMat* reconObs );
	
	// release memories;
	void Release();

	
private:
	//core structures	
    int nInitialObs;	          // the number of samples for initial observation
	int nKDim;	                  //target dimension
	float ff;                       // forgotten factor
	float m_updateThreshold;
	int probUpdate;

	// principle axes, each axes is a column of the matrix
	CvMat *pPrinAxes;
	// to rotate previous points into updated coordinates
	CvMat* pUnitBasis;
	//mean of input histogram
	CvMat * pMean;
	// a storage of inpute color histogrrams for initlization
	CvMat* pObservation;
	// updating pca needing
	CvMat* pVeeHeader;// prinaxes v's header
	CvMat* pInterCol;//intermediate coloumn		
	CvMat * pYou;

    // aux parameters
	__int64 nSampleCounter;// to record the samples received

};
