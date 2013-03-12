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
                const DImage &flow, const DImage &rflow,
                int start=2, int end=10, double na=15);
    int cluster(DImage &centers, DImage &layers,
                const DImage &features,
                int width, int height, int start=2, int end=10, double na=15);

    static double mydist(double *p1, double *p2, int start, int end);
private:
    const static int sWidth = 2;// spatial info width
    const static int iWidth = 1;// image info width
    const static int fWidth = 2;// flow info width
};

#endif
