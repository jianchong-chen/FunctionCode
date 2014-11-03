#include "stdafx.h"

#include "LHParas.h"
#include <math.h>
#include "LHCCIPCA.h"
#include <iostream>
#include <fstream>
#define WINDOW_SIZE 1000000
using namespace std;

void CLHCCIPCA::Initial(CvMat *pFirstObs, int targetDim, int batchFrame, float forgetFattor)
{
	if(pFirstObs->cols < targetDim) return;
	// see header files for more comment 
	nInitialObs = batchFrame;// the number of samples for initial observation
	nKDim = targetDim;//target dimension
	ff = forgetFattor;
	m_updateThreshold = 0;

	// each colomn is an principle axes
	pPrinAxes = cvCreateMat(nKDim,pFirstObs->cols,pFirstObs->type);
	pUnitBasis = cvCreateMat(nKDim,pFirstObs->cols,pFirstObs->type);	
	pMean = cvCreateMat(1,pFirstObs->cols,pFirstObs->type);		
	pObservation = cvCreateMat(nInitialObs,pFirstObs->cols,pFirstObs->type);
	pVeeHeader = cvCreateMatHeader(1,pFirstObs->cols,pFirstObs->type);// prinaxes v's header
	pInterCol = cvCreateMat(1,pFirstObs->cols,pFirstObs->type);//intermediate coloumn		
	pYou=cvCreateMat(1,pFirstObs->cols,pFirstObs->type);
	
	nSampleCounter = 0L; 

	return;
}
void CLHCCIPCA::Initial2(void)
{
	return;
}
void CLHCCIPCA::Release()
{
	cvReleaseMat(&pPrinAxes);
	cvReleaseMat(&pMean);
	cvReleaseMat(&pYou);

	cvReleaseMatHeader(&pVeeHeader);
	cvReleaseMat(&pInterCol);
	cvReleaseMat(&pUnitBasis);
}

bool CLHCCIPCA::UpdatePCA( const CvMat* pCurrentObs, float &updateErr )
{ 
	assert(pCurrentObs->cols == pPrinAxes->cols);
	assert(pCurrentObs->type == pPrinAxes->type);

	//m_updateThreshold = updateThreshold;

	__int64 t1, t2;
	float tF = cvGetTickFrequency(); 

	if(nSampleCounter < nInitialObs)
	{
		cvGetRow(pObservation, pVeeHeader, nSampleCounter);
		cvCopy(pCurrentObs, pVeeHeader);

		nSampleCounter++;
		return false;
	}
	else if(nSampleCounter == nInitialObs)
	{
		// run pca
		int minValue = (pObservation->rows - pObservation->cols) > 0 ? pObservation->cols : pObservation->rows;
		CvMat* pTempEigenVect=cvCreateMat(minValue,pCurrentObs->cols,pCurrentObs->type);
		CvMat* pEigenValue=cvCreateMat(minValue,1,pCurrentObs->type);

		cvCalcPCA(pObservation,pMean,pEigenValue,pTempEigenVect,CV_PCA_DATA_AS_ROW);
		// because open cv output each principle vector per row, so we transfer and use first nKDim of them.
		for(int j=0;j<nKDim;j++)
		{
			for(int i=0;i<pCurrentObs->rows;i++)
			{
				CV_MAT_ELEM(*pPrinAxes,float,j,i) = CV_MAT_ELEM(*pTempEigenVect,float,j,i);
			}
		}

		double eigenValue = 0;
		cvCopy(pPrinAxes, pUnitBasis);
		//cvTranspose(pUnitBasis,pUnitBasisT);
		//// norm to unit basis		
		for(int i=0;i<pPrinAxes->rows;i++)
		{
			cvGetRow(pPrinAxes,pVeeHeader,i);			
			eigenValue = CV_MAT_ELEM(*pEigenValue,float,i,0);
			if (_isnan(eigenValue) || eigenValue < 0.00000001)
			{
				break;
			}
			cvScale(pVeeHeader,pVeeHeader,eigenValue);
		}
		
		cvReleaseMat(&pObservation);
		cvReleaseMat(&pEigenValue);
		cvReleaseMat(&pTempEigenVect);

		nSampleCounter++;
		return false;
	}
	else
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		//   follows are for normal state
		////////////////////////////////////////////////////////////////////////////////////////////////////////

		float tF = cvGetTickFrequency();
	
		//__int64 N = nSampleCounter;		
		//float alpha=(1.0f+ff)/N;
		float alpha = LH_LEARN_RATE;  // original is 0.001
		
		float normYou = 0;
		float normVee = 0;
		float uV = 0;

		nSampleCounter++;
		if(nSampleCounter>WINDOW_SIZE)
		{
			nSampleCounter = nInitialObs+1;
		}

		cvSub(pCurrentObs,pMean,pYou);		
		if (_isnan(CV_MAT_ELEM(*pYou, float, 0, 0)))
		{
			return false;
		}

		for(int i=0;i<nKDim;i++)
		{
			// update \v_i+1
			cvGetRow(pPrinAxes,pVeeHeader,i);
			normVee = (float)cvNorm(pVeeHeader);
			//float small_value = pow(10,-6);
			if (_isnan(normVee) || normVee < 0.000001)
			{
				normVee = 0.000001;
				//break;
			}

			uV = (float)cvDotProduct(pYou,pVeeHeader);
			if (_isnan(uV) || abs(uV) < 0.000001)
			{
				uV = normVee;
				//break;
			}
			cvScale(pYou,pInterCol,uV);	
			cvAddWeighted(pVeeHeader, 1.0f-alpha, pInterCol, alpha/normVee, 0, pVeeHeader);

			// update \u_i+1
			uV = (float)cvDotProduct(pYou,pVeeHeader);			
			normVee = (float)cvNorm(pVeeHeader);
			if (_isnan(normVee) || _isnan(uV))
			{
				break;
			}
			if(normVee < 0.00001)
			{
				normVee = 0.00001;
			}
			uV = -1*uV/(normVee*normVee);
			cvAddWeighted(pYou, 1.0f, pVeeHeader, uV, 0, pYou);
		}
		// update mean
		updateErr = cvNorm(pYou);
		if (_isnan(updateErr))
		{
			return false;
		}
		cvAddWeighted(pMean, 1-alpha, pCurrentObs, alpha, 0, pMean);

		//// norm to unit basis	
		cvCopy(pPrinAxes,pUnitBasis);
		for(int i=0;i<pPrinAxes->rows;i++)
		{
			cvGetRow(pUnitBasis,pVeeHeader,i);
			normVee = (float)cvNorm(pVeeHeader);
			if (_isnan(normVee) || normVee < 0.000001)
			{
				normVee = 0.000001;
				//break;
			}
			cvScale(pVeeHeader,pVeeHeader,1/normVee);
		}

		return true;
	}
}

int CLHCCIPCA::MatchSubBrick(const CvMat* pCurrentObs, float& errRecon, float m_Threshold_Divid)
{
	__int64 t10, t2;
	float tF = cvGetTickFrequency();
	float uV = 0;
	float normVee = 0;
	float errRecon_Rate = 0;
	float meanNorm = 0;
	CvScalar sumYou;
	
	assert(pCurrentObs->cols == pPrinAxes->cols);
	assert(pCurrentObs->type == pPrinAxes->type);

	t10 = cvGetTickCount();
	cvSub(pCurrentObs, pMean, pYou);
	normVee = cvNorm(pYou);
	if (_isnan(normVee))
	{
		return true;
	}

	// Computing the matching Threshold that is adaptive to the intensity of the average brick
	sumYou = cvSum(pYou);
	meanNorm = cvNorm(pMean);

	float Divid_Num = LH_DIVID_NUM; 
	
	if(sumYou.val[0] > 480)
	{
		m_updateThreshold = meanNorm/Divid_Num + sumYou.val[0]/LH_DIVID_NUM_1;
	}
	else if(sumYou.val[0] < -560)
	{
		m_updateThreshold = meanNorm/Divid_Num + sumYou.val[0]/LH_DIVID_NUM_2;
	}
	else
	{
		m_updateThreshold = meanNorm/Divid_Num;// + sumYou.val[0]/100;
	}
	if(m_updateThreshold < LH_MIN_THRESHOLD)
	{
		m_updateThreshold = LH_MIN_THRESHOLD;  //16
	}

	// Computing the Reconstructed Error of the current brick against the background model
	for(int i=0; i<nKDim; i++)
	{
		cvGetRow(pUnitBasis,pVeeHeader,i);
		uV = (float)cvDotProduct(pYou,pVeeHeader);	
		if (_isnan(uV))
		{
			break;
		}
		cvAddWeighted(pYou, 1.0f, pVeeHeader, -uV, 0, pYou);
	}
	errRecon = cvNorm(pYou);
	if (_isnan(errRecon))
	{
		return true;
	}
	errRecon_Rate = errRecon/normVee;

	t2 = cvGetTickCount();
	//cout<<"In_Match(): "<<(t2-t10)/(tF*1000)<<"ms"<<endl;
	
	if( errRecon > m_updateThreshold /*|| errRecon_Rate < 0.49*/ )   // mThreshold default value is 80
	{
		return 0;
	}
	else if( errRecon < 4 )  // 重建误差太小，“太是背景”，也不更新
	{
		return 2;
	}
	else
	{
		return 1;
	}
}


bool CLHCCIPCA::ReconBrick( const CvMat* pCurrentObs, CvMat* reconObs )
{
	// add the mean
	cvAddWeighted(pMean, 1, reconObs, 1, 0, reconObs);
	return true;
}