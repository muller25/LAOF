#ifndef _MotionLayers_H
#define _MotionLayers_H

#include "Image.h"
#include "ML.h"
#include "Maths.h"
#include "ImageProcess.h"

class MotionLayers
{
public:
    void spatialInfo(DImage &info, const DImage &im);
    void imInfo(DImage &info, const DImage &im);
    void flowInfo(DImage &info, const DImage &flow);
    int cluster(DImage &centers, DImage &layers, const DImage &im,
                const DImage &flow, const DImage &rflow);
    
    template <class T>
    static double mydist(T *p1, T *p2, int start, int end);
private:
    const static int sWidth = 2;// spatial info width
    const static int iWidth = 1;// image info width
    const static int fWidth = 2;// flow info width
};

template <class T>
double MotionLayers::mydist(T *p1, T *p2, int start, int end)
{
    double dist = 0;
    int idx = 0;
    
    // spatial info
    // dist += dist2(p1, p2, idx, idx+sWidth);
    // idx += sWidth;
    
    // image info
    // dist += dist2(p1, p2, idx, idx+iWidth);
    // idx += iWidth;
    
    // flow info
    dist += dist2(p1, p2, idx, idx+fWidth);
    dist += (1 - similarity(p1, p2, idx, idx+fWidth));
    idx += fWidth;
    
    // reverse flow info
    dist += dist2(p1, p2, idx, idx+fWidth);
    dist += (1 - similarity(p1, p2, idx, idx+fWidth));

    return dist;
}

#endif
