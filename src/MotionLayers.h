#ifndef _MotionLayers_H
#define _MotionLayers_H

#include <cv.h>

#include "Image.h"

typedef std::Mat_<double> DMat;

class MotionLayers
{
public:
    void flow2feature(DMat &features, const DImage &im, const DImage &flow);
    int clusters(UCImage &labels, const DMat &features);
};

#endif
