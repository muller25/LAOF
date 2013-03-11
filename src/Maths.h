#ifndef _Maths_H
#define _Maths_H

#include <cmath>
#include <vector>

#include "Image.h"

#ifndef ESP
#define ESP 1e-6
#endif

const double PI = atan(1.0) * 4;

// res = m1 * m2 * m3
template <class T, class T1, class T2, class T3>
inline void multiply(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2, const Image<T3> &m3)
{
    assert(m1.match3D(m2) && m2.match3D(m3));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = (T)((double)m1[i] * (double)m2[i] * (double)m3[i]);
}

// res = m1 * m2
template <class T, class T1, class T2>
inline void multiply(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = (T)((double)m1[i] * (double)m2[i]);
}

// res = m1 * val
template <class T, class T1, class T2>
inline void multiply(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = (T)((double)m1[i] * val);
}

// res *= m1
template <class T, class T1>
inline void multiply(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] *= m1[i];
}

// res *= val
template <class T, class T1>
inline void multiply(Image<T> &res, const T1 val)
{
    assert(res.ptr() != NULL);

    for (int i = 0; i < res.nElements(); ++i)
        res[i] *= val;
}

// res = m1 / m2
template <class T, class T1, class T2>
inline void divide(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] / m2[i];
}

// res = m1 / val
template <class T, class T1, class T2>
inline void divide(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    assert(val != 0);

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] / val;
}

// res /= m1
template <class T, class T1>
inline void divide(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));

    for (int i = 0; i < res.nElements(); ++i)
        res[i] /= m1[i];
}

// res /= val
template <class T, class T1>
inline void divide(Image<T> &res, const T1 val)
{
    assert(val != 0 && res.ptr() != NULL);

    for (int i = 0; i < res.nElements(); ++i)
        res[i] /= val;
}

// res = m1 + m2
template <class T, class T1, class T2>
inline void add(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] + m2[i];
}

// res = m1 + val
template <class T, class T1, class T2>
inline void add(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] + val;
}

// res += m1
template <class T, class T1>
inline void add(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));

    for (int i = 0; i < res.nElements(); ++i)
        res[i] += m1[i];
}

// res += val
template <class T, class T1>
inline void add(Image<T> &res, const T1 val)
{
    assert(res.ptr() != NULL);
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] += val;
}

// res = m1 - m2
template <class T, class T1, class T2>
inline void substract(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] - m2[i];
}

// res = m1 - val
template <class T, class T1, class T2>
inline void substract(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] - val;
}

// res -= m1
template <class T, class T1>
inline void substract(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));

    for (int i = 0; i < res.nElements(); ++i)
        res[i] += m1[i];
}

// res -= val
template <class T, class T1>
inline void substract(Image<T> &res, const T1 val)
{
    assert(res.ptr() != NULL);
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] -= val;
}

// res -= val * m1
template <class T, class T1>
inline void substract(Image<T> &res, const Image<T1> &m1, const double val)
{
    assert(res.ptr() != NULL);
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] -= m1[i] * val;
}

// res = m1 * val + m2 * val + val3
template <class T, class T1, class T2>
inline void addWeighted(Image<T> &res, const Image<T1> &m1, double val1,
                        const Image<T2> &m2, double val2, T val3=0)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < res.nElements(); ++i)
        res[i] = m1[i] * val1 + m2[i] * val2 + val3;
}

template <class T>
inline T enforceRange(const T &num, const T &minVal, const T &maxVal)
{
    return std::min(std::max(num, minVal), maxVal);
}

template<class T, class T1>
inline void biInterpolate(T *res, const Image<T1> &m, const double x, const double y)
{
    int width = m.nWidth(), height = m.nHeight(), channels = m.nChannels();
    int nx, ny, offset;
    double s, dx, dy;
    int ix = x;
    int iy = y;
    T1 *pm = m.ptr();

    dx = std::max(std::min(x-ix, 1.), 0.);
    dy = std::max(std::min(y-iy, 1.), 0.);
    for (int hh = 0; hh <= 1; hh++)
    {
        ny = enforceRange(iy+hh, 0, height-1);
        for (int ww = 0; ww <= 1; ww++)
        {
            nx = enforceRange(ix+ww, 0, width-1);
            s = fabs(1-hh-dy) * fabs(1-ww-dx);
            offset = (ny * width + nx) * channels;

            for (int k = 0; k < channels; ++k)
                res[k] += pm[offset+k] * s;
        }
    }
}

template <class T, class T1>
void hfiltering(Image<T> &res, const Image<T1> &src, double *filter, int fsize)
{
    int width = src.nWidth(), height = src.nHeight(), channels = src.nChannels();
    T1 *ps = src.ptr();
    
    res.create(width, height, channels);
    T *pr = res.ptr();
    double factor;
    int roffset, soffset, ww;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            roffset = (h * width + w) * channels;
            for (int f = -fsize; f <= fsize; ++f)
            {
                ww = enforceRange(f+w, 0, width-1);
                factor = filter[f+fsize];
                soffset = (h * width + ww) * channels;
                for (int k = 0; k < channels; ++k)
                    pr[roffset+k] += ps[soffset+k] * factor;
            }
        }
    }
}

template <class T, class T1>
void vfiltering(Image<T> &res, const Image<T1> &src, double *filter, int fsize)
{
    int width = src.nWidth(), height = src.nHeight(), channels = src.nChannels();
    T1 *ps = src.ptr();
    
    res.create(width, height, channels);
    T *pr = res.ptr();
    double factor;
    int roffset, soffset, hh;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            roffset = (h * width + w) * channels;
            for (int f = -fsize; f <= fsize; ++f)
            {
                hh = enforceRange(f+h, 0, height-1);
                factor = filter[f+fsize];
                soffset = (hh * width + w) * channels;

                for (int k = 0; k < channels; ++k)
                    pr[roffset+k] += (double)ps[soffset+k] * factor;
            }
        }
    }
}

template <class T, class T1>
void filtering(Image<T> &res, const Image<T1> &src, double *hfilter, int hfsize,
               double *vfilter, int vfsize)
{
    Image<T> tmp;
    hfiltering(tmp, src, hfilter, hfsize);
    vfiltering(res, tmp, vfilter, vfsize);
}

// central difference
template <class T, class T1>
void gradX(Image<T> &res, const Image<T1> &src)
{
    static double filter[] = {1./12, -8./12, 0, 8./12, -1./12};
    hfiltering(res, src, filter, 2);
}

// central difference
template <class T, class T1>
void gradY(Image<T> &res, const Image<T1> &src)
{
    static double filter[] = {1./12, -8./12, 0, 8./12, -1./12};
    vfiltering(res, src, filter, 2);
}

template <class T, class T1>
void grad1st(Image<T> &gx, Image<T> &gy, const Image<T1> &src)
{
    static double filter[] = {1./12, -8./12, 0, 8./12, -1./12};
    hfiltering(gx, src, filter, 2);
    vfiltering(gy, src, filter, 2);
}

template <class T, class T1>
void gradXX(Image<T> &res, const Image<T1> &src)
{
    static double filter[] = {-1./12, 16./12, -30./12, 16./12, -1./12};
    hfiltering(res, src, filter, 2);
}

template <class T, class T1>
void gradYY(Image<T> &res, Image<T1> &src)
{
    static double filter[] = {-1./12, 16./12, -30./12, 16./12, -1./12};
    vfiltering(res, src, filter, 2);
}

template <class T, class T1>
void gradXY(Image<T> &res, const Image<T1> &src)
{
    // double filter[3][3] = {
    //     {.25, 0, -.25},
    //     {0,   0,    0},
    //     {-.25,0,  .25}
    // };

    Image<T> tmp;
    tmp = gradX(tmp, src);
    gradY(res, tmp);
}

// forward difference for u and backward difference for phi'
template<class T, class T1, class T2>
void weighted_lap(Image<T> &lap, const Image<T1> &flow, const Image<T2> &weight)
{
    assert(flow.match3D(weight) && flow.nChannels() == 1);
    
    int width = flow.nWidth(), height = flow.nHeight();
    T1 *pf = flow.ptr();
    T2 *pw = weight.ptr();

    lap.create(width, height);

    T *pl = lap.ptr();
    int offset;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            if (w < width-1)
                pl[offset] += pw[offset] * (pf[offset+1]-pf[offset]);
            if (w > 0)
                pl[offset] -= pw[offset-1] * (pf[offset]-pf[offset-1]);
            if (h < height-1)
                pl[offset] += pw[offset] * (pf[offset+width]-pf[offset]);
            if (h > 0)
                pl[offset] -= pw[offset-width] * (pf[offset]-pf[offset-width]);
        }
    }
}

// laplacian operator for spatial-temporal
// pf = previous flow, f(t-1), cf = current flow, f(t), nf = next flow, f(t+1)
template<class T, class F>
void weighted_lap3(Image<T> &lap, const Image<F> &pf, const Image<F> &cf,
                   const Image<F> &nf, const Image<T> &pw, const Image<T> &cw)
{
    assert(pf.match3D(cf) && cf.match3D(nf) && nf.match3D(pw) && pw.match3D(cw) &&
           cw.nChannels() == 1);

    weighted_lap(lap, cf, cw);
    for (int i = 0; i < lap.nElements(); ++i)
        lap[i] += cw[i] * (nf[i]-cf[i]) - pw[i] * (cf[i]-pf[i]);
}

template <class T, class T1>
void collapse(Image<T> &dst, const Image<T1> &src)
{
    int width = src.nWidth(), height = src.nHeight(), channels = src.nChannels();

    if (channels == 1)
    {
        src.copyTo(dst);
        return;
    }

    T1 *ps = src.ptr();
    
    dst.create(width, height);
    T *pd = dst.ptr();
    double tmp;
    int doffset, soffset;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            doffset = h * width + w;
            soffset = doffset * channels;
            tmp = 0;
            for (int k = 0; k < channels; ++k)
                tmp += ps[soffset  + k];
            
            pd[doffset] = (T)((double)tmp / channels);
        }
    }
}

template <class T, class T1>
inline bool equal(const Image<T> &m1, const Image<T1> &m2)
{
    assert(m1.match3D(m2));
    for (int i = 0; i < m1.nElements(); ++i)
        if (fabs(m1[i] - m2[i]) > ESP) return false;

    return true;
}

template <class T>
void split(std::vector< Image<T> > &arr, Image<T> &m)
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

// merge channels
template <class T>
void mergec(Image<T> &m, std::vector< Image<T> > &arr)
{
    int vsize = arr.size();
    int nchannels = arr[0].nChannels();
    for (int i = 1; i < vsize; ++i)
    {
        assert(arr[0].match2D(arr[i]));
        nchannels += arr[i].nChannels();
    }

    int size = arr[0].nSize(), channels, o, mo, c = 0;

    m.create(arr[0].nWidth(), arr[0].nHeight(), nchannels);
    for (int i = 0; i < vsize; ++i)
    {
        channels = arr[i].nChannels();
        for (int s = 0; s < size; ++s)
        {
            mo = s * nchannels;
            o = s * channels;
            for (int k = 0; k < channels; ++k)
                m[mo+k+c] = arr[i][o+k];
        }
        c += channels;
    }
}

template <class T>
void mergew(Image<T> &m, std::vector< Image<T> > &arr)
{
    int height = arr[0].nHeight(), channels = arr[0].nChannels();
    int nwidth = arr[0].nWidth(), ww, width, o, mo;
    int vsize = arr.size();
    for (int v = 1; v < vsize; ++v)
    {
        assert(height == arr[v].nHeight() && channels == arr[v].nChannels());
        nwidth += arr[v].nWidth();
    }

    m.create(nwidth, height, channels);
    ww = 0;
    for (int v = 0; v < vsize; ++v)
    {
        width = arr[v].nWidth();
        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; ++w)
            {
                o = (h * width + w) * channels;
                mo = (h * nwidth + w + ww) * channels;
                for (int k = 0; k < channels; ++k)
                    m[mo+k] = arr[v][o+k];
            }
        }
        
        ww += arr[v].nWidth();
    }
}

template <class T>
void mergeh(Image<T> &m, std::vector< Image<T> > &arr)
{
    int width = arr[0].nWidth(), channels = arr[0].nChannels();
    int nheight = arr[0].nHeight(), hh, height, o, mo;
    int vsize = arr.size();
    for (int v = 1; v < vsize; ++v)
    {
        assert(width == arr[v].nHeight() && channels == arr[v].nChannels());
        nheight += arr[v].nHeight();
    }

    m.create(width, nheight, channels);
    hh = 0;
    for (int v = 0; v < vsize; ++v)
    {
        height = arr[v].nHeight();
        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; ++w)
            {
                o = (h * width + w) * channels;
                mo = ((h + hh) * width + w) * channels;
                for (int k = 0; k < channels; ++k)
                    m[mo+k] = arr[v][o+k];
            }
        }
        
        hh += arr[v].nHeight();
    }
}

template <class T>
double averError(Image<T> &m1, Image<T> &m2)
{
    assert(m1.match3D(m2));
    
    double error = 0;
    for (int i = 0; i < m1.nElements(); ++i)
        error += fabs(m1[i] - m2[i]);

    error /= m1.nElements();
    return error;
}

// normalize each channels seperately
template <class T>
void normalizeChannels(Image<T> &res, const Image<T> &m)
{
    int width = m.nWidth(), height = m.nHeight(), channels = m.nChannels();
    int size = width*height, offset;
    T *mins = new T[channels], *maxs = new T[channels];

    res.create(width, height, channels);
    // init
    for (int k = 0; k < channels; ++k)
    {
        mins[k] = m[k];
        maxs[k] = m[k];
    }

    // get max and min of each channel
    for (int i = 1; i < size; ++i)
    {
        offset = i * channels;
        for (int k = 0; k < channels; ++k)
        {
            mins[k] = std::min(m[offset+k], mins[k]);
            maxs[k] = std::max(m[offset+k], maxs[k]);
        }
    }

    // normalize
    for (int i = 0; i < size; ++i)
    {
        offset = i * channels;
        for (int k = 0; k < channels; ++k)
            res[offset+k] = (m[offset+k] - mins[k]) / (maxs[k] - mins[k]);
    }
    
    delete []mins;
    delete []maxs;
}

// normalize
template <class T>
void normalize(Image<T> &res, const Image<T> &m)
{
    T maxm = m.max(), minm = m.min();
    
    res.create(m.nWidth(), m.nHeight(), m.nChannels());
    for (int i = 0; i < m.nElements(); ++i)
        res[i] = (m[i] - minm) / (maxm - minm);
}

// Euclidean distance
// calculate range [start, end)
template <class T>
inline double dist2(T *p1, T *p2, int start, int end)
{
    double d = 0;
    for (int i = start; i < end; ++i)
        d += (p1[i] - p2[i]) * (p1[i] - p2[i]);

    return sqrt(d);
}

// city block distance
template <class T>
inline double dist1(T *p1, T *p2, int start, int end)
{
    double dist = 0;
    for (int i = start; i < end; ++i)
        dist += fabs(p1[i] - p2[i]);

    return dist;
}

// Minkowski distance
template <class T>
inline double distp(T *p1, T *p2, double p, int start, int end)
{
    double dist = 0;

    for (int i = start; i < end; ++i)
        dist += pow(fabs(p1[i]-p2[i]), p);
    
    return pow(dist, 1./p);
}

// vector similarity from start to end, end is exclusive
// using cos as measurement
template <class T>
double similarity(T *v1, T *v2, int start, int end)
{
    T a, b, ab;
    a = b = ab = 0;
    for (int i = start; i < end; ++i)
    {
        a += v1[i] * v1[i];
        b += v2[i] * v2[i];
        ab += v1[i] * v2[i];
    }

    a = sqrt(a), b = sqrt(b);
    if (fabs(a) < ESP) // a == 0
        return b;
    
    if (fabs(b) < ESP) // b == 0
        return a;

    return (double)ab / (a * b);
}

// O(1) time rectangular sum, independent of window size
// use border replicate for pixel out of boundary
// each element is mean of window 2*wsize+1, 2*hsize+1
template <class T, class T1>
void rectSum(Image<T> &dst, const Image<T1> &src, int wsize=1, int hsize=1)
{
    int width = src.nWidth(), height = src.nHeight(), channels = src.nChannels();
    int offset, idx, last, wc = width * channels;
    Image<T> cumSum(width, height, channels), tmp(width, height, channels);
    
    dst.create(width, height, channels);
    
    // cumulative sum over y axis
    for (int w = 0; w < width; ++w)
    {
        idx= w * channels;
        for (int k = 0; k < channels; ++k)
            cumSum[idx+k] = (hsize+1) * src[idx+k];
    }

    for (int h = 1; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            idx = (h * width + w) * channels;
            for (int k = 0; k < channels; ++k)
                cumSum[idx+k] = src[idx+k] + cumSum[idx-wc+k];
        }
    }
    
    offset = hsize * width * channels; // corresponding row in cumSum
    last = (height-1) * width;         // last row
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            idx = (h * width + w) * channels;
            for (int k = 0; k < channels; ++k)
            {
                if (h <= hsize)
                    tmp[idx+k] = cumSum[idx+offset+k] - src[w*channels+k]*h;
                else if (h > hsize && h + hsize < height)
                    tmp[idx+k] = cumSum[idx+offset+k] - cumSum[idx-offset-wc+k];
                else // h + hsize >= height
                    tmp[idx+k] = cumSum[(last+w)*channels+k] - cumSum[idx-offset-wc+k] +
                        (h+hsize-height+1) * src[(last+w)*channels+k];
            }
        }
    }

    // cumulative sum over x axis
    for (int h = 0; h < height; ++h)
    {
        idx = h * width * channels;
        for (int k = 0; k < channels; ++k)
            cumSum[idx+k] = (wsize+1) * tmp[idx+k];
    }

    for (int h = 0; h < height; ++h)
    {
        for (int w = 1; w < width; ++w)
        {
            idx = (h * width + w) * channels;
            for (int k = 0; k < channels; ++k)
                cumSum[idx+k] = tmp[idx+k] + cumSum[idx-channels+k];
        }
    }
    
    offset = wsize * channels;
    last = (width-1) * channels; // last column
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            idx = (h * width + w) * channels;
            for (int k = 0; k < channels; ++k)
            {
                if (w <= wsize)
                    dst[idx+k] = cumSum[idx+offset+k] - tmp[h*width*channels+k]*w;
                else if (w > wsize && w + wsize < width)
                    dst[idx+k] = cumSum[idx+offset+k] - cumSum[idx-offset-channels+k];
                else // w + wsize >= width
                    dst[idx+k] = cumSum[h*width*channels+last+k] - cumSum[idx-offset-channels+k] +
                        (w+wsize-width+1) * tmp[(h*width*channels)+last+k];
            }
        }
    }
}

// brute force to calculate rectangular sum
template <class T, class T1>
void rectSumBF(Image<T> &dst, const Image<T1> &src, int wsize=1, int hsize=1)
{
    int width = src.nWidth(), height = src.nHeight(), channels = src.nChannels();
    int hh, ww, offset;
    
    dst.create(width, height, channels);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            for (int hs = -hsize; hs <= hsize; ++hs)
            {
                hh = enforceRange(h+hs, 0, height-1);
                for (int ws = -wsize; ws <= wsize; ++ws)
                {
                    ww = enforceRange(w+ws, 0, width-1);
                    for (int k = 0; k < channels; ++k)
                        dst[offset+k] += src[(hh*width+ww)*channels+k];
                }
            }
        }
    }
}

template <class T>
void randFill(Image<T> &m, T minV, T maxV)
{
    assert(m.ptr() != NULL);
    srand(time(NULL));
    
    T range = maxV - minV;
    for (int i = 0; i < m.nElements(); ++i)
    {
        if (m.isFloat())
            m[i] = minV + (double)rand() / RAND_MAX * range;
        else
            m[i] = minV + rand() % (int)(range+1);
    }
}

#endif
