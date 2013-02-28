#include "GaussianPyramid.h"

void GaussianPyramid::ConstructPyramid(const Mat &im, double ratio, int minWidth)
{
	// the ratio cannot be arbitrary numbers
	if(ratio > 0.98 || ratio < 0.4)
		ratio = 0.75;
    
	// first decide how many levels
	levels = log((double)minWidth / im.cols) / log(ratio);

	if (pPyr != NULL)
		delete []pPyr;
    
	pPyr = new Mat[levels];
    im.copyTo(pPyr[0]);

	int n = log(0.25) / log(ratio); 
	double baseSigma = (1 / ratio - 1), nSigma = baseSigma * n, sigma, rate;
    Mat tmp;
    
	for(int i = 1; i < levels; i++)
	{
		if(i <= n)
		{
			sigma = baseSigma * i;
            GaussianBlur(im, tmp, Size(0, 0), sigma, sigma, BORDER_REPLICATE);
            resize(tmp, pPyr[i], Size(0, 0), pow(ratio, i), pow(ratio, i));
		}
		else
		{
            GaussianBlur(pPyr[i-n], tmp, Size(0, 0), nSigma, nSigma, BORDER_REPLICATE);
            rate = pow(ratio, i) * (double)im.cols / tmp.cols;
            resize(tmp, pPyr[i], Size(0, 0), rate, rate);
		}
	}
}

void GaussianPyramid::ConstructPyramid(const Mat &im, int minWidth)
{
    levels = ceil(log((double)minWidth / im.cols) / log(0.5));

    if (pPyr != NULL)
        delete []pPyr;

    pPyr = new Mat[levels];
    im.copyTo(pPyr[0]);

    for (int i = 1; i < levels; i++)
        pyrDown(pPyr[i-1], pPyr[i], Size(0, 0), BORDER_REPLICATE);
}
