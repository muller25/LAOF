#ifndef _HIST_H
#define _HIST_H

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
}

// oeraitation * strength
template <class T, class T1, class T2>
void calcOMHist(Image<T> &hist, const Image<T1> &m, const Image<T2> &mask,
                int bins=18)
{
    assert(m.nChannels() == 2 && mask.nChannels() == 1);
    double interval = 360. / bins, theta;
    T1 x, y;
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
    }
}

#endif
