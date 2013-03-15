#ifndef _LAOF_H
#define _LAOF_H

#include "Image.h"
#include <vector>

class LAOF
{
public:
    static void EM(std::vector<DImage> &u, std::vector<DImage> &v,
                   std::vector<DImage> &ur, std::vector<DImage> &vr,
                   std::vector<DImage> &masks, const std::vector<DImage> &im,
                   int curIdx, int nIters=3);

    static void transferLabels(DImage &dst, const DImage &src, const DImage &flow);
};

#endif
