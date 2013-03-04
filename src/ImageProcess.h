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
void gaussianBlur(Image<T> &dst, const Image<T1> &src, double sigma, int fsize)
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

#endif
