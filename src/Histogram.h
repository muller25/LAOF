#ifndef _Histogram_H
#define _Histogram_H

#include "Image.h"

template <class T, class T1, class T2>
void calcHist(Image<T> &hist, const Image<T1> &m, const Image<T2> &mask,
              T1 minV, T1 maxV, int bins=16)
{
    assert(maxV > minV && mask.nChannels() == 1);
    
    int channels = m.nChannels(), size = m.nSize();
    double interval = (double)(maxV - minV + 1) / bins;
    int idx;

    hist.create(bins * channels, 1);
    for (int i = 0; i < size; ++i)
    {
        // not in mask
        if (mask.isZero(i)) continue;
        
        for (int k = 0; k < channels; ++k)
        {
            idx = (m[i*channels+k] - minV) / interval;
            hist[k*bins+idx]++;
        }
    }

    // normalize
    T total;
    for (int k = 0; k < channels; ++k)
    {
        total = 0;
        for (int b = 0; b < bins; ++b)
            total += hist[k*bins+b];

        for (int b = 0; b < bins; ++b)
            hist[k*bins+b] /= total;
    }
}

template <class T, class T1>
T getHistVal(Image<T> &hist, const Image<T1> &m, int idx, T1 minV, T1 maxV)
{
    int channels = m.nChannels();
    int bins = hist.nWidth() / channels;
    double interval = (double)(maxV - minV + 1) / bins;

    T val = 0;
    int imVal;
    for (int k = 0; k < channels; ++k)
    {
        imVal = (m[idx*channels+k] - minV) / interval;
        val += hist[k*bins+imVal];
    }

    return val;
}

// oeraitation * strength
template <class T, class T1, class T2>
void calcOMHist(Image<T> &hist, const Image<T1> &m, const Image<T2> &mask,
                int bins=18)
{
    assert(m.nChannels() == 2 && mask.nChannels() == 1);
    double interval = 360. / bins, theta;
    T1 x, y;
    T total = 0;
    int idx;
    
    hist.create(bins, 1);
    for (int i = 0; i < m.nSize(); ++i)
    {
        if (mask.isZero(i)) continue;
   
        x = m[i*2];
        y = m[i*2+1];

        // arc to degree
        theta = (atan2(y, x) + PI) / PI * 180.;
        idx = theta / interval;
        hist[idx] += theta * sqrt(x*x + y*y);
        total += hist[idx];
    }

    for (int i = 0; i < bins; ++i)
        hist[idx] /= total;
}

template <class T, class T1>
T getOMHistVal(Image<T> &hist, const Image<T1> &m, int idx)
{
    assert(m.nChannels() == 2);
    
    int bins = hist.nWidth();
    double interval = 360./bins;
    T1 x = m[idx*2];
    T1 y = m[idx*2+1];
    double theta = (atan2(y, x) + PI) / PI * 180.;
    int pos = theta / interval;
    return hist[pos];
}

#endif
