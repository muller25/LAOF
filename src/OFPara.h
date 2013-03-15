#ifndef _OFPara_H
#define _OFPara_H

#include "Image.h"
#include <vector>

class OFPara
{
public:
    DImage im1, im2, mask1, mask2;
    DImage u1, v1, u2, v2;
    int lux, luy, rbx, rby;

    OFPara(){lux = luy = rbx = rby = -1;}
    virtual ~OFPara()
    {
        im1.release();
        im2.release();
        mask1.release();
        mask2.release();
        u1.release();
        v1.release();
        u2.release();
        v2.release();
        lux = luy = rbx = rby = -1;
    }
    
    static void split(std::vector<OFPara> &paras, int labels,
                      const DImage &im1, const DImage &im2,
                      const DImage &l1, const DImage &l2);

    static void restore(DImage &flow, DImage &rflow, const std::vector<OFPara> &paras,
                        int width, int height);

    static void restore(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                        const std::vector<OFPara> &paras, int width, int height);
private:
    const static int border = 20;
};

#endif
