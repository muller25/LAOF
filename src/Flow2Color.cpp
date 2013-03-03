#include "Flow2Color.h"

void flow2color(UCImage &flowImg, UCImage &idxImg, const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const double thresh = 10e9;
    
    int width = u.nWidth(), height = u.nHeight(), offset, ioffset;
    double *pu = u.ptr(), *pv = v.ptr();

    // fix unknown flow
    idxImg.create(width, height);
    uchar *pIdx = idxImg.ptr();
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            if (pu[offset] >= thresh || pu[offset] <= -thresh)
            {
                pu[offset] = 0;
                pIdx[ioffset] = 255;
            }
            
            if (pv[offset] >= thresh || pv[offset] <= -thresh)
            {
                pv[offset] = 0;
                pIdx[ioffset] = 255;
            }
        }
    }
    
    // find min and max
    DImage rad(width, height);
    double *pr = rad.ptr();
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            pr[offset] = sqrt(pu[offset]*pu[offset] + pv[offset]*pv[offset]);
        }
    }
    double maxrad = rad.max();
    
    printf("max flow: %.4f flow range: ", maxrad);
    printf("u = %.3f .. %.3f; ", rad.max(), u.min(),  u.max());
    printf("v = %.3f .. %.3f\n", v.min(), v.max());

    DImage nu, nv;
    if (fabs(maxrad) < ESP)
    {
        u.copyTo(nu);
        v.copyTo(nv);
    }
    else
    {
        divide(nu, u, maxrad);
        divide(nv, v, maxrad);
    }
    
    computeColor(flowImg, nu, nv);
}

void computeColor(UCImage &im, const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const double PI = atan(1) * 4;

    int width = u.nWidth(), height = u.nHeight();
    double *pu = u.ptr(), *pv = v.ptr();

    DImage rad(width, height);
    double *pr = rad.ptr();

    IImage k0(width, height), k1(width, height);
    int *pk = k0.ptr(), *pk1 = k1.ptr();

    DImage f(width, height);    
    double *pf = f.ptr(), fk;

    DImage wheel = makeColorWheel();
    double *pw = wheel.ptr();
    int ncols = wheel.nHeight(), offset;

    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;

            // rad = sqrt(u^2 + v^2)
            pr[offset] = pv[offset]*pv[offset] + pu[offset]*pu[offset];
            pr[offset] = sqrt(pr[offset]);
            
            // a = atan2(-v, -u) / pi
            fk = atan2(-pv[offset], -pu[offset]) / PI;

            // -1~1 mapped to 1~ncols
            fk = (fk+1) / 2 * (ncols-1) + 1;

            // 1, 2, ..., ncols
            pk[offset] = (int)floor(fk);

            if (pk[offset] == ncols)
                pk1[offset] = 1;
            else
                pk1[offset] = pk[offset] + 1;

            pf[offset] = fk - pk[offset];
        }
    }

    im.create(width, height, ncols);
    uchar *pi = im.ptr();
    double col0, col1, col;
    bool idx;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            idx = (pr[offset] <= 1);
        
            for (int k = 0; k < ncols; ++k)
            {
                col0 = pw[pk[offset] * width + k] / 255;
                col1 = pw[pk1[offset] * width + k] / 255;
                col = (1 - pf[offset]) * col0 + pf[offset] * col1;

                // increase saturation with radius
                if (idx) col = 1 - pr[offset] * (1 - col);
                // out of range
                else col = col * 0.75;

                pi[offset*ncols + k] = (uchar)floor(255 * col);
            }
        }
    }
}

DImage makeColorWheel()
{
    const int RY = 15;
    const int YG = 6;
    const int GC = 4;
    const int CB = 11;
    const int BM = 13;
    const int MR = 6;
    const int width = 3;
    int height = RY + YG + GC + CB + BM + MR;
    
    DImage wheel(width, height);
    double *pw = wheel.ptr();

    // RY
    int start = 0;
    for (int h = start; h < start + RY; ++h)
    {
        pw[h * width] = 255;
        pw[h * width + 1] = floor(255 * h / RY);
    }
    start += RY;
    
    // YG
    for (int h = start; h < start + YG; ++h)
    {
        pw[h * width] = 255 - floor(255 * (h-start) / YG);
        pw[h * width + 1] = 255;
    }
    start += YG;

    // GC
    for (int h = start; h < start + GC; ++h)
    {
        pw[h * width + 1] = 255;
        pw[h * width + 2] = floor(255 * (h-start) / GC);
    }
    start += GC;

    // CB
    for (int h = start; h < start + CB; ++h)
    {
        pw[h * width + 1] = 255 - floor(255 * (h-start) / CB);
        pw[h * width + 2] = 255;
    }
    start += CB;

    // BM
    for (int h = start; h < start + BM; ++h)
    {
        pw[h * width + 2] = 255;
        pw[h * width] = floor(255 * (h-BM) / BM);
    }
    start += BM;
    
    // MR
    for (int h = start; h < start + MR; ++h)
    {
        pw[h * width + 2] = 255 - floor(255 * (h-MR) / MR);
        pw[h * width] = 255;
    }

    return wheel;
}
