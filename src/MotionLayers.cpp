#include "MotionLayers.h"

void MotionLayers::flow2feature(DMat &features, const DImage &im, const DImage &flow)
{
    const int imBins = 36;
    const int omBins = 16;
    assert(im.match2D(flow));
    int size = im.nSize();
    
    features.create(size, imBins+flowBins);
    
}

int MotionLayers::clusters(UCImage &labels, const DMat &features)
{
}
