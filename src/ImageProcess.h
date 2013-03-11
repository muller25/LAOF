#ifndef _ImageProcess_H
#define _ImageProcess_H

#include "Image.h"
#include "Maths.h"

template <class T, class T1>
void imresize(Image<T> &dst, const Image<T1> &src, double factor)
{
    assert(src.ptr() != NULL && factor > 0);

    int srcWidth = src.nWidth();
    int srcHeight = src.nHeight();
    int channels = src.nChannels();
    int dstWidth = (double)srcWidth * factor;
    int dstHeight = (double)srcHeight * factor;
    double x, y;
    
    dst.create(dstWidth, dstHeight, channels);
    T *ptr = dst.ptr();

    for(int h = 0; h < dstHeight; ++h)
    {
        for (int w = 0; w < dstWidth; ++w)
		{
			y = (double)(h+1) / factor - 1;
			x = (double)(w+1) / factor - 1;

			// bilinear interpolation
			biInterpolate(ptr + (h * dstWidth + w) * channels, src, x, y);
		}
    }
}

template <class T, class T1>
void imresize(Image<T> &dst, const Image<T1> &src, int dstWidth, int dstHeight)
{
    assert(src.ptr() != NULL && dstWidth > 0 && dstHeight > 0);
    
    int srcWidth = src.nWidth(), srcHeight = src.nHeight(), channels = src.nChannels();
    dst.create(dstWidth, dstHeight, channels);    

	double xRatio = (double)dstWidth / srcWidth;
	double yRatio = (double)dstHeight / srcHeight;
	double x,y;

    T *ptr = dst.ptr();
	for (int h = 0; h < dstHeight; ++h)
    {
		for(int w = 0; w < dstWidth; ++w)
		{
			x = (double)(w+1) / xRatio - 1;
			y = (double)(h+1) / yRatio - 1;

			// bilinear interpolation
			biInterpolate(ptr + (h * dstWidth + w) * channels, src, x, y);
		}
    }
}

template <class T, class T1>
void desuarate(Image<T> &res, const Image<T1> &im)
{
    assert (im.ptr() != NULL);

    int width = im.nWidth(), height = im.nHeight(), channels = im.nChannels();
    T1 *pi = im.ptr();
    
    res.create(width, height);
    T *pr = res.ptr();
    int go, of;
    
    // desuarate BGR -> Gray
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            go = h * width + w;
            of = go * channels;
            pr[go] = (double)pi[of+2]*.299 + (double)pi[of+1]*.587 + (double)pi[of]*.114;
        }
    }
}

template <class T, class T1>
void im2double(Image<T> &dst, const Image<T1> &src)
{
    assert(dst.isFloat());
    
    if (src.isFloat())
    {
        src.copyTo(dst);
        return;
    }

    dst.create(src.nWidth(), src.nHeight(), src.nChannels());
    for (int i = 0; i < src.nElements(); i++)
        dst[i] = (T) src[i] / 255;
}

// filters
// gaussian smoothing
template <class T, class T1>
void gSmooth(Image<T> &dst, const Image<T1> &src, double sigma, int fsize)
{
    // constructing the 1D gaussian filter
	double* gFilter = new double[fsize * 2 + 1];
    assert(gFilter != NULL);
    
	double sum = 0;
	sigma = sigma * sigma * 2;

	for(int i = -fsize; i <= fsize; i++)
	{
		gFilter[i + fsize] = exp(-(double)(i * i) / sigma);
		sum += gFilter[i + fsize];
	}

	for(int i = 0; i < 2*fsize+1; i++)
		gFilter[i] /= sum;
    
	// apply filtering
	filtering(dst, src, gFilter, fsize, gFilter, fsize);

	delete []gFilter;
}

// O(1) time box filter
template <class T, class T1>
void BoxFilter(Image<T> &dst, const Image<T1> &src, int wsize=3, int hsize=3)
{
    rectSum(dst, src, wsize, hsize);

    double factor = 1. / ((2 * wsize + 1) * (2 * hsize + 1));
    multiply(dst, factor);
}

// O(1) time guided image filter, for single channel images only
// r is for local window radius, try 2, 4, 8
// eps is for regularization, try 0.1^2, 0.2^2, 0.4^2
template <class T, class T1, class T2>
void GuidedFilter(Image<T> &dst, const Image<T1> &im, const Image<T2> &guide,
                  int r=4, double eps=0.04)
{
    assert(guide.match3D(im) && im.nChannels() == 1);

    DImage meanI, meanp, meanIp, covIp, tmp, meanI2, varI, a, b, meanA, meanB;

    BoxFilter(meanI, guide, r, r); // meanI = mean (guide)
    BoxFilter(meanp, im, r, r);    // meanp = mean (im)

    multiply(tmp, guide, im);
    BoxFilter(meanIp, tmp, r, r); // meanIp = mean (guide .* im)
    
    // this is the covariance of (I, p) in each local patch.
    multiply(tmp, meanI, meanp);
    substract(covIp, meanIp, tmp); // covIp = meanIp - meanI * meanp
    
    multiply(tmp, guide, guide);
    BoxFilter(meanI2, tmp, r, r); // meanI2 = mean (guide * guide)

    multiply(tmp, meanI, meanI);  
    substract(varI, meanI2, tmp); // varI = meanI2 - meanI * meanI
    add(varI, eps);

    divide(a, covIp, varI);       // a = covIp ./ (varI + eps)
    multiply(meanI, a);
    substract(b, meanp, meanI);   // b = meanp - a .* meanI
    
    BoxFilter(meanA, a, r, r);    // a = mean (a)
    BoxFilter(meanB, b, r, r);    // b = mean (b)

    multiply(dst, meanA, guide);  
    add(dst, meanB);              // dst = meanA * guide + meanB
}

// O(1) time guided image filter
// for multi-channel guide image and single channel input image
template <class T, class T1, class T2>
void GuidedFilterColor(Image<T> &dst, const Image<T1> &im, const Image<T2> &guide,
                       int r=4, double eps=0.04)
{
    assert(im.match2D(guide) && guide.nChannels() == 3 && im.nChannels() == 1);

    DImage meanIr, meanIg, meanIb, meanp, meanIpr, meanIpg, meanIpb, tmp;
    DImage covIpr, covIpg, covIpb;
    std::vector<DImage> g3;

    split(g3, guide); // BGR
    
    BoxFilter(meanIb, g3[0], r, r); // meanI = mean (guide)
    BoxFilter(meanIg, g3[1], r, r);
    BoxFilter(meanIr, g3[2], r, r);
    BoxFilter(meanp, im, r, r);    // meanp = mean (im)
    
    multiply(tmp, g3[0], im);
    BoxFilter(meanIpb, tmp, r, r); // meanIpb = mean (g3[b])

    multiply(tmp, g3[1], im);
    BoxFilter(meanIpg, tmp, r, r); // meanIpg = mean (g3[g])

    multiply(tmp, g3[2], im);
    BoxFilter(meanIpr, tmp, r, r); // meanIpr = mean (g3[r])
    
    // covariance of (I, p) in each local patch.
    multiply(tmp, meanIr, meanp);
    substract(covIpr, meanIpr, tmp); // covIpr = meanIpr - meanIr * meanp

    multiply(tmp, meanIg, meanp);
    substract(covIpg, meanIpg, tmp); // covIpg = meanIpg - meanIg * meanp

    multiply(tmp, meanIb, meanp);
    substract(covIpb, meanIpb, tmp); // covIpb = meanIpb - meanIb * meanp

    DImage varIr2, varIrg, varIrb, varIg2, varIgb, varIb2;

    // varIr2 = mean (g[r] .* g[r]) - meanIr .* meanIr
    multiply(tmp, g3[2], g3[2]);
    BoxFilter(varIr2, tmp, r, r);
    multiply(tmp, meanIr, meanIr);
    substract(varIr2, tmp); 

    // varIrg = mean (g[r] .* g[g]) - meanIr .* meanIg
    multiply(tmp, g3[2], g3[1]);
    BoxFilter(varIrg, tmp, r, r);
    multiply(tmp, meanIr, meanIg);
    substract(varIrg, tmp); 

    // varIrb = mean (g[r] .* g[b]) - meanIr .* meanIb
    multiply(tmp, g3[2], g3[0]);
    BoxFilter(varIrb, tmp, r, r);
    multiply(tmp, meanIr, meanIb);
    substract(varIrb, tmp); 

    // varIg2 = mean (g[g] .* g[g]) - meanIg .* meanIg
    multiply(tmp, g3[1], g3[1]);
    BoxFilter(varIg2, tmp, r, r);
    multiply(tmp, meanIg, meanIg);
    substract(varIg2, tmp); 

    // varIgb = mean (g[g] .* g[b]) - meanIg .* meanIb
    multiply(tmp, g3[1], g3[0]);
    BoxFilter(varIgb, tmp, r, r);
    multiply(tmp, meanIg, meanIb);
    substract(varIgb, tmp); 

    // varIb2 = mean (g[b] .* g[b]) - meanIb .* meanIb
    multiply(tmp, g3[0], g3[0]);
    BoxFilter(varIb2, tmp, r, r);
    multiply(tmp, meanIb, meanIb);
    substract(varIb2, tmp);
    
    // variance of I in each local patch: the matrix Sigma in Eqn (14).
    // Note the variance in each local patch is a 3x3 symmetric matrix:
    //           rr, rg, rb
    //   Sigma = rg, gg, gb
    //           rb, gb, bb
    int width = im.nWidth(), height = im.nHeight(), offset;
    DImage sigma(3, 3), covIp(3, 1), a(width, height, 3), b, inv;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            sigma[0] = varIr2[offset] + eps;
            sigma[1] = varIrg[offset];
            sigma[2] = varIrb[offset];
            sigma[3] = varIrg[offset];
            sigma[4] = varIg2[offset] + eps;
            sigma[5] = varIgb[offset];
            sigma[6] = varIrb[offset];
            sigma[7] = varIgb[offset];
            sigma[8] = varIb2[offset] + eps;
            sigma.invTo(inv);

            // covIp = [covIpr, covIpg, covIpb]
            covIp[0] = covIpr[offset];
            covIp[1] = covIpg[offset];
            covIp[2] = covIpb[offset];

            // a = covIp * inv(sigma + eps * eyes(3))
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    a[offset+j] += covIp[k] * inv[k*3+j];
        }
    }

    std::vector<DImage> a3;
    split(a3, a);

    // b = meanp - a[0] * meanIr - a[1] * meanIg - a[2] * meanIb
    multiply(tmp, a3[0], meanIr);
    substract(b, meanp, tmp);

    multiply(tmp, a3[1], meanIg);
    substract(b, tmp);

    multiply(tmp, a3[2], meanIb);
    substract(b, tmp); 

    // dst = mean(a[0]) * g[r] + mean(a[1]) * g[b] + mean(a[2]) * g[b] + mean(b)
    BoxFilter(tmp, a3[0], r, r);
    multiply(dst, tmp, g3[2]);

    BoxFilter(tmp, a3[1], r, r);
    multiply(tmp, g3[1]);
    add(dst, tmp);

    BoxFilter(tmp, a3[2], r, r);
    multiply(tmp, g3[0]);
    add(dst, tmp);

    BoxFilter(tmp, b, r, r);
    add(dst, tmp);
}

#endif
