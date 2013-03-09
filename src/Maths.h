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
void multiply(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2, const Image<T3> &m3)
{
    assert(m1.match3D(m2) && m2.match3D(m3));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = (T)((double)m1[i] * (double)m2[i] * (double)m3[i]);
}

// res = m1 * m2
template <class T, class T1, class T2>
void multiply(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = (T)((double)m1[i] * (double)m2[i]);
}

// res = m1 * val
template <class T, class T1, class T2>
void multiply(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = (T)((double)m1[i] * val);
}

// res *= m1
template <class T, class T1>
void multiply(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] *= m1[i];
}

// res *= val
template <class T, class T1>
void multiply(Image<T> &res, const T1 val)
{
    assert(res.ptr() != NULL);

    for (int i = 0; i < res.nElements(); ++i)
        res[i] *= val;
}

// res = m1 / m2
template <class T, class T1, class T2>
void divide(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] / m2[i];
}

// res = m1 / val
template <class T, class T1, class T2>
void divide(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    assert(val != 0);

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] / val;
}

// res /= m1
template <class T, class T1>
void divide(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));

    for (int i = 0; i < res.nElements(); ++i)
        res[i] /= m1[i];
}

// res /= val
template <class T, class T1>
void divide(Image<T> &res, const T1 val)
{
    assert(val != 0 && res.ptr() != NULL);

    for (int i = 0; i < res.nElements(); ++i)
        res[i] /= val;
}

// res = m1 + m2
template <class T, class T1, class T2>
void add(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] + m2[i];
}

// res = m1 + val
template <class T, class T1, class T2>
void add(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] + val;
}

// res += m1
template <class T, class T1>
void add(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));

    for (int i = 0; i < res.nElements(); ++i)
        res[i] += m1[i];
}

// res += val
template <class T, class T1>
void add(Image<T> &res, const T1 val)
{
    assert(res.ptr() != NULL);
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] += val;
}

// res = m1 - m2
template <class T, class T1, class T2>
void substract(Image<T> &res, const Image<T1> &m1, const Image<T2> &m2)
{
    assert(m1.match3D(m2));

    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] - m2[i];
}

// res = m1 - val
template <class T, class T1, class T2>
void substract(Image<T> &res, const Image<T1> &m1, const T2 val)
{
    res.create(m1.nWidth(), m1.nHeight(), m1.nChannels());
    for (int i = 0; i < m1.nElements(); ++i)
        res[i] = m1[i] - val;
}

// res -= m1
template <class T, class T1>
void substract(Image<T> &res, const Image<T1> &m1)
{
    assert(res.match3D(m1));

    for (int i = 0; i < res.nElements(); ++i)
        res[i] += m1[i];
}

// res -= val
template <class T, class T1>
void substract(Image<T> &res, const T1 val)
{
    assert(res.ptr() != NULL);
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] -= val;
}

// res -= val * m1
template <class T, class T1>
void substract(Image<T> &res, const Image<T1> &m1, const double val)
{
    assert(res.ptr() != NULL);
    
    for (int i = 0; i < res.nElements(); ++i)
        res[i] -= m1[i] * val;
}

// res = m1 * val + m2 * val + val3
template <class T, class T1, class T2>
void addWeighted(Image<T> &res, const Image<T1> &m1, double val1,
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
bool equal(Image<T> &m1, Image<T1> &m2)
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

template <class T>
void merge(Image<T> &m, std::vector< Image<T> > &arr)
{
    assert(!arr.empty());
    int width = arr[0].nWidth(), height = arr[0].nHeight(), channels = arr.size(), offset;

    m.create(width, height, channels);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; h < width; ++w)
        {
            offset = h * width + w;
            for (int k = 0; k < channels; ++k)
                m[offset*channels+k] = arr[k][offset];
        }
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

// Euclidean distance
template <class T>
inline double dist2(T *p1, T *p2, int size)
{
    double d = 0;
    for (int i = 0; i < size; ++i)
        d += (p1[i] - p2[i]) * (p1[i] - p2[i]);

    return sqrt(d);
}

#endif
