#ifndef _ImageProcess_H
#define _ImageProcess_H

#include "Image.h"
#include "Maths.h"

// cover labels onto image
template <class I, class L>
void coverLabels(Image<I> &res, const Image<I> &im, const Image<L> &labels)
{
    assert(im.match2D(labels));

    const double factor = 0.3;
    int width = im.nWidth(), height = im.nHeight(), channels = im.nChannels();
    int lchannels = labels.nChannels();
    
    res.create(width, height, channels);
    if (lchannels == 1)
    {
        int ioffset;
        for(int i = 0; i < labels.nSize(); ++i)
        {
            ioffset = i * channels;
            for (int k = 0; k < channels; ++k)
            {
                res[ioffset+k] = im[ioffset+k]*factor;
            
                if (typeid(I) == typeid(L))
                    res[ioffset+k] += labels[i]*(1-factor);
                else if (im.isFloat() && !labels.isFloat())
                    res[ioffset+k] += (double)labels[i]/255.*(1-factor);
                else
                    res[ioffset+k] += labels[i]*(1-factor)*255;
            }
        }
    } else {
        assert(lchannels == channels);
        
        for (int i = 0; i < labels.nElements(); ++i)
        {
            res[i] = im[i]*factor;
            
            if (typeid(I) == typeid(L))
                res[i] += labels[i]*(1-factor);
            else if (im.isFloat() && !labels.isFloat())
                res[i] += (double)labels[i]/255.*(1-factor);
            else
                res[i] += labels[i]*(1-factor)*255;
        }
    }
}

template <class I, class F>
void warpImage(Image<I> &warp, const Image<I> &im1, const Image<I> &im2,
               const Image<F> &u, const Image<F> &v)
{
    assert(im1.match3D(im2) && u.match3D(v) && im1.match2D(u));
    
    int width = im1.nWidth(), height = im1.nHeight(), channels = im1.nChannels();
    double nx, ny;
    int offset;

    warp.create(width, height, channels);    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            nx = w + u[offset];
            ny = h + v[offset];
            offset *= channels;
            
            if (nx < 0 || nx > width-1 || ny < 0 || ny > height-1)
            {
                for (int k = 0; k < channels; k++)
                    warp[offset+k] = im1[offset+k];

                continue;
            }

            biInterpolate(warp.ptr()+offset, im2, nx, ny);
        }
    }
}

template <class I, class F>
void warpImage(Image<I> &warp, const Image<I> &im1, const Image<I> &im2,
               const Image<F> &flow)
{
    assert(im1.match3D(im2) && im1.match2D(flow) && flow.nChannels() == 2);
    
    int width = im1.nWidth(), height = im1.nHeight(), channels = im1.nChannels();
    double nx, ny;
    int offset;

    warp.create(width, height, channels);    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            nx = w + flow[offset*2];
            ny = h + flow[offset*2+1];
            offset *= channels;
            
            if (nx < 0 || nx > width-1 || ny < 0 || ny > height-1)
            {
                for (int k = 0; k < channels; k++)
                    warp[offset+k] = im1[offset+k];

                continue;
            }

            biInterpolate(warp.ptr()+offset, im2, nx, ny);
        }
    }
}

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

template <class T>
void split(std::vector< Image<T> > &arr, const Image<T> &m)
{
    assert(m.ptr() != NULL);

    int width = m.nWidth(), height = m.nHeight(), channels = m.nChannels(), offset;
    
    arr.clear();
    for (int k = 0; k < channels; ++k)
        arr.push_back(Image<T>(width, height));
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            for (int k = 0; k < channels; ++k)
                arr[k][offset] = m[offset*channels+k];
        }
    }
}

// 显示mask覆盖的图像区域。mask的取值应该只有0，1，且显示mask=1的区域
template <class I, class M>
void cut(Image<I> &res, const Image<I> &im, const Image<M> &mask)
{
    assert(im.match2D(mask) && mask.nChannels() == 1);

    int offset, channels = im.nChannels();
    res.create(im.nWidth(), im.nHeight(), channels);
    for (int i = 0; i < im.nSize(); ++i)
    {
        if (mask.isZero(i)) continue;

        offset = i * channels;
        for (int k = 0; k < channels; ++k)
            res[offset + k] = im[offset + k];
    }
}

// 显示矩形覆盖的图像区域。lux, luy为矩形左上角坐标，rbx,rby为矩形右下角坐标
template <class T>
void cut(Image<T> &dst, const Image<T> &src, int lux, int luy, int rbx, int rby)
{
    assert(lux <= rbx && luy <= rby && src.ptr() != NULL);

    int width = rbx - lux + 1;
    int height = rby - luy + 1;
    int sWidth = src.nWidth(), channels = src.nChannels();
    int sOf, dOf;
    dst.create(width, height, channels);

    for (int h = luy; h <= rby; ++h)
    {
        for (int w = lux; w <= rbx; ++w)
        {
            sOf = (h * sWidth + w) * channels;
            dOf = ((h-luy) * width + (w-lux)) * channels;
            for (int k = 0; k < channels; ++k)
                dst[dOf+k] = src[sOf+k];
        }
    }
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
    std::vector< Image<T2> > g3;

    split(g3, guide); // BGR
    BoxFilter(meanIb, g3[B], r, r); // meanI = mean (guide)
    BoxFilter(meanIg, g3[G], r, r);
    BoxFilter(meanIr, g3[R], r, r);
    BoxFilter(meanp, im, r, r);    // meanp = mean (im)
    
    multiply(tmp, g3[B], im);
    BoxFilter(meanIpb, tmp, r, r); // meanIpb = mean (g3[b] * im)

    multiply(tmp, g3[G], im);
    BoxFilter(meanIpg, tmp, r, r); // meanIpg = mean (g3[g] * im)

    multiply(tmp, g3[R], im);
    BoxFilter(meanIpr, tmp, r, r); // meanIpr = mean (g3[r] * im)
    
    // covariance of (I, p) in each local patch.
    multiply(tmp, meanIr, meanp);
    substract(covIpr, meanIpr, tmp); // covIpr = meanIpr - meanIr * meanp

    multiply(tmp, meanIg, meanp);
    substract(covIpg, meanIpg, tmp); // covIpg = meanIpg - meanIg * meanp

    multiply(tmp, meanIb, meanp);
    substract(covIpb, meanIpb, tmp); // covIpb = meanIpb - meanIb * meanp

    DImage varIr2, varIrg, varIrb, varIg2, varIgb, varIb2;

    // varIr2 = mean (g[r] .* g[r]) - meanIr .* meanIr
    multiply(tmp, g3[R], g3[R]);
    BoxFilter(varIr2, tmp, r, r);
    multiply(tmp, meanIr, meanIr);
    substract(varIr2, tmp); 

    // varIrg = mean (g[r] .* g[g]) - meanIr .* meanIg
    multiply(tmp, g3[R], g3[G]);
    BoxFilter(varIrg, tmp, r, r);
    multiply(tmp, meanIr, meanIg);
    substract(varIrg, tmp); 

    // varIrb = mean (g[r] .* g[b]) - meanIr .* meanIb
    multiply(tmp, g3[R], g3[B]);
    BoxFilter(varIrb, tmp, r, r);
    multiply(tmp, meanIr, meanIb);
    substract(varIrb, tmp); 

    // varIg2 = mean (g[g] .* g[g]) - meanIg .* meanIg
    multiply(tmp, g3[G], g3[G]);
    BoxFilter(varIg2, tmp, r, r);
    multiply(tmp, meanIg, meanIg);
    substract(varIg2, tmp); 

    // varIgb = mean (g[g] .* g[b]) - meanIg .* meanIb
    multiply(tmp, g3[G], g3[B]);
    BoxFilter(varIgb, tmp, r, r);
    multiply(tmp, meanIg, meanIb);
    substract(varIgb, tmp); 

    // varIb2 = mean (g[b] .* g[b]) - meanIb .* meanIb
    multiply(tmp, g3[B], g3[B]);
    BoxFilter(varIb2, tmp, r, r);
    multiply(tmp, meanIb, meanIb);
    substract(varIb2, tmp);
    
    // variance of I in each local patch: the matrix Sigma in Eqn (14).
    // Note the variance in each local patch is a 3x3 symmetric matrix:
    //           rr, rg, rb
    //   Sigma = rg, gg, gb
    //           rb, gb, bb
    int width = im.nWidth(), height = im.nHeight(), offset;
    DImage sigma(3, 3), a(width, height, 3), b, inv;
    
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
            // a = covIp * inv(sigma + eps * eyes(3))
            for (int k = 0; k < 3; ++k)
                a[offset*3+k] = covIpr[offset] * inv[k] +
                                covIpg[offset] * inv[3+k] +
                                covIpb[offset] * inv[6+k];
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

    BoxFilter(dst, b, r, r);

    // dst = mean(a[0]) * g[r] + mean(a[1]) * g[g] + mean(a[2]) * g[b] + mean(b)
    BoxFilter(tmp, a3[0], r, r);
    multiply(tmp, g3[R]);
    add(dst, tmp);
    
    BoxFilter(tmp, a3[1], r, r);
    multiply(tmp, g3[G]);
    add(dst, tmp);

    BoxFilter(tmp, a3[2], r, r);
    multiply(tmp, g3[B]);
    add(dst, tmp);
}

#endif
