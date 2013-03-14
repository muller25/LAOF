#ifndef _LAOF_H
#define _LAOF_H

#include "OpticalFlow.h"
#include "MotionLayers.h"
#include "Image.h"

#include "ImageIO.h"
#include "Flow2Color.h"

class LAOF
{
public:
    static void EM(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                   DImage &mask1, DImage &mask2,
                   const DImage &im1, const DImage &im2,
                   int nIters=3);

    static void createCenterByLabel(DImage &centers, int clusters,
                                    const DImage &labels, const DImage &samples);

    static void transferLabel(DImage &dst, const DImage &src, const DImage &flow);
};

#endif
