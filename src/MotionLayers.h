#ifndef _MotionLayers_H
#define _MotionLayers_H

#include "Image.h"
#include "ML.h"
#include "Maths.h"

class MotionLayers
{
public:
    void im2feature(DImage &features, const DImage &im);
    void flow2feature(DImage &features, const DImage &flow, const DImage &rflow);
    int flowCluster(DImage &centers, DImage &layers, const DImage &im,
                    const DImage &flow, const DImage &rflow);
    
    template <class T>
    static double mdist(T *p1, T *p2, int size);
};

template <class T>
double MotionLayers::mdist(T *p1, T *p2, int size)
{
    double dist = dist2(p1, p2, size);

    // vector similarity
    T a, b, ab;
    a = b = ab = 0;
    for (int i = 0; i < size; ++i)
    {
        a += p1[i] * p1[i];
        b += p2[i] * p2[i];
        ab += p1[i] * p2[i];
    }

    a = sqrt(a), b = sqrt(b);
    if (fabs(a) < ESP) // a == 0
        dist += b;
    else if (fabs(b) < ESP) // b == 0
        dist += a;
    else
        dist += (1 - ab / (a * b));
    
    return dist;
}

#endif
