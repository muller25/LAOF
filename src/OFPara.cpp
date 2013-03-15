#include "OFPara.h"
#include "ImageProcess.h"
#include "Maths.h"

void OFPara::split(std::vector<OFPara> &paras, int labels,
                   const DImage &im1, const DImage &im2,
                   const DImage &l1, const DImage &l2)
{
    assert(im1.match3D(im2) && l1.match3D(l2) && l1.nChannels() == 1);

    paras.clear();

    int width = im1.nWidth(), height = im1.nHeight();
    int l1_lux, l1_luy, l1_rbx, l1_rby, l2_lux, l2_luy, l2_rbx, l2_rby;
    int lux, luy, rbx, rby, offset;
    OFPara p;
    DImage tmp1, tmp2;

    for (int l = 0; l < labels; ++l)
    {
        l1_lux = width, l1_luy = height, l1_rbx = -1, l1_rby = -1;
        l2_lux = width, l2_luy = height, l2_rbx = -1, l2_rby = -1;
        tmp1.create(width, height);
        tmp2.create(width, height);

        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; ++w)
            {
                offset = h * width + w;
                if (fabs(l1[offset] - l) < 0.5)
                {
                    tmp1[offset] = 1;
                    l1_lux = std::min(l1_lux, w);
                    l1_luy = std::min(l1_luy, h);
                    l1_rbx = std::max(l1_rbx, w);
                    l1_rby = std::max(l1_rby, h);
                }

                if (fabs(l2[offset] - l) < 0.5)
                {
                    tmp2[offset] = 1;
                    l2_lux = std::min(l2_lux, w);
                    l2_luy = std::min(l2_luy, h);
                    l2_rbx = std::max(l2_rbx, w);
                    l2_rby = std::max(l2_rby, h);
                }
            }
        }

        lux = enforceRange(std::min(l1_lux, l2_lux) - border, 0, width-1);
        luy = enforceRange(std::min(l1_luy, l2_luy) - border, 0, height-1);
        rbx = enforceRange(std::max(l1_rbx, l2_rbx) + border, 0, width-1);
        rby = enforceRange(std::max(l1_rby, l2_rby) + border, 0, height-1);

        cutRect(p.mask1, tmp1, lux, luy, rbx, rby);
        cutRect(p.mask2, tmp2, lux, luy, rbx, rby);
        cutRect(p.im1, im1, lux, luy, rbx, rby);
        cutRect(p.im2, im2, lux, luy, rbx, rby);
        p.lux = lux, p.luy = luy, p.rbx = rbx, p.rby = rby;
        paras.push_back(p);
    }
}

void OFPara::restore(DImage &flow, DImage &rflow, const std::vector<OFPara> &paras,
                     int width, int height)
{
    assert(!paras.empty());
    
    flow.create(width, height, 2);
    rflow.create(width, height, 2);

    int offset, suboff, subWidth, subHeight;
    for (size_t l = 0; l < paras.size(); ++l)
    {
        subWidth = paras[l].rbx - paras[l].lux + 1;
        subHeight = paras[l].rby - paras[l].luy + 1;
        for (int hh = 0; hh < subHeight; ++hh)
        {
            for (int ww = 0; ww < subWidth; ++ww)
            {
                offset = (hh + paras[l].luy) * width + ww + paras[l].lux;
                suboff = hh * subWidth + ww;
                offset *= 2;

                if (fabs(paras[l].mask1[suboff] - 1) < ESP)
                {
                    flow[offset] = paras[l].u1[suboff];
                    flow[offset+1] = paras[l].v1[suboff];
                }
                if (fabs(paras[l].mask2[suboff] - 1) < ESP)
                {
                    rflow[offset] = paras[l].u2[suboff];
                    rflow[offset+1] = paras[l].v2[suboff];
                }
            }
        }
    }
}

void OFPara::restore(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                     const std::vector<OFPara> &paras,
                     int width, int height)
{
    assert(!paras.empty());
    
    u1.create(width, height);
    v1.create(width, height);
    u2.create(width, height);
    v2.create(width, height);

    int offset, suboff, subWidth, subHeight;
    for (size_t l = 0; l < paras.size(); ++l)
    {
        subWidth = paras[l].rbx - paras[l].lux + 1;
        subHeight = paras[l].rby - paras[l].luy + 1;
        for (int hh = 0; hh < subHeight; ++hh)
        {
            for (int ww = 0; ww < subWidth; ++ww)
            {
                offset = (hh + paras[l].luy) * width + ww + paras[l].lux;
                suboff = hh * subWidth + ww;

                if (fabs(paras[l].mask1[suboff] - 1) < ESP)
                {
                    u1[offset] = paras[l].u1[suboff];
                    v1[offset] = paras[l].v1[suboff];
                }

                if (fabs(paras[l].mask2[suboff] - 1) < ESP)
                {
                    u2[offset] = paras[l].u2[suboff];
                    v2[offset] = paras[l].v2[suboff];
                }
            }
        }
    }
}

