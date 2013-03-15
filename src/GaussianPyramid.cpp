#include "GaussianPyramid.h"

void GaussianPyramid::build(const DImage &im, double ratio, int minWidth, bool smooth)
{
	// the ratio cannot be arbitrary numbers
	if(ratio > 0.98 || ratio < 0.4)
		ratio = 0.75;
    
	// first decide how many levels
	levels = log((double)minWidth / im.nWidth()) / log(ratio);
    assert(levels > 0);
    
	if (pPyr != NULL)
		delete []pPyr;
    
	pPyr = new DImage[levels];
    assert(pPyr != NULL);
    im.copyTo(pPyr[0]);

	int n = log(0.25) / log(ratio); 
	double baseSigma = (1/ratio - 1);
    double nSigma = baseSigma * n;
    double sigma, rate;
    DImage tmp;
    
	for (int i = 1; i < levels; ++i)
	{
		if(i <= n)
		{
			sigma = baseSigma * i;
            if (smooth)
            {
                gSmooth(tmp, im, sigma, sigma*3);
                imresize(pPyr[i], tmp, pow(ratio, i));
            }
            else imresize(pPyr[i], im, pow(ratio, i));
		}
		else
		{
            if (smooth)
            {
                gSmooth(tmp, pPyr[i-n], nSigma, nSigma*3);
                rate = pow(ratio, i) * (double)im.nWidth() / tmp.nWidth();
                imresize(pPyr[i], tmp, rate);
            } else {
                rate = pow(ratio, i) * (double)im.nWidth() / pPyr[i-n].nWidth();
                imresize(pPyr[i], pPyr[i-n], rate);
            }
		}
	}
}
