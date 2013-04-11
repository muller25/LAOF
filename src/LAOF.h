#ifndef _LAOF_H
#define _LAOF_H

#include "Image.h"
#include <vector>

class LAOF
{
public:
    static void EM(std::vector<DImage> &u, std::vector<DImage> &v,
                   std::vector<DImage> &ur, std::vector<DImage> &vr,
                   std::vector<DImage> &mask, const std::vector<DImage> &im,
                   int cur, int nIters=3);

    static void segment(DImage &centers1, DImage &centers2,
                        DImage &layers1, DImage &layers2, int &nlabels,
                        const DImage &im1, const DImage &im2,
                        const DImage &u, const DImage &v,
                        const DImage &ur, const DImage &vr);
    
    static void transferLabels(DImage &dst, const DImage &src, const DImage &flow);
    static void intFlow(DImage &dst, const DImage &src);
    static char out[256];
};

#endif
