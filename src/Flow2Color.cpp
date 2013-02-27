#include "Flow2Color.h"

void flow2color(Mat &u, Mat &v, Mat &flowImg, Mat &idxImg)
{
    assert(matchAll(u, v) && u.type() == CV_64F);

    const double thresh = 10e9;
    
    int rows = u.rows, cols = u.cols;
    double *pu = (double *)u.data, *pv = (double *)v.data;
    int step = u.step / sizeof(double);
    int r, c, offset;

    // fix unknown flow
    Mat unknownIdx = Mat::zeros(rows, cols, CV_64F);
    double *pIdx = (double *)unknownIdx.data;
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c;
            if (pu[offset] >= thresh || pu[offset] <= -thresh)
            {
                pu[offset] = 0;
                pIdx[offset] = 1;
            }
            
            if (pv[offset] >= thresh || pv[offset] <= -thresh)
            {
                pv[offset] = 0;
                pIdx[offset] = 1;
            }
        }
    }
    unknownIdx.convertTo(idxImg, CV_8U);
    
    // find min and max
    double minu, maxu, minv, maxv, maxrad;
    Mat rad = u * u + v * v;
    minMaxIdx(u, &minu, &maxu);
    minMaxIdx(v, &minv, &maxv);
    minMaxIdx(rad, NULL, &maxrad);
    maxrad = sqrt(maxrad);

    printf("max flow: %.4f flow range: u = %.3f .. %.3f;", maxrad, minu, maxu);
    printf("v = %.3f .. %.3f\n", minv, maxv);

    Mat nu, nv;
    if (fabs(maxrad) < ESP)
    {
        nu = u;
        nv = v;
    }
    else
    {
        nu = u / maxrad;
        nv = v / maxrad;
    }
    
    flowImg = computeColor(nu, nv);
}

Mat computeColor(Mat &u, Mat &v)
{
    assert(matchAll(u, v) && u.type() == CV_64F);

    const double PI = atan(1) * 4;        

    int rows = u.rows, cols = u.cols;
    int step = u.step / sizeof(double), offset;
    double *pu = (double *)u.data, *pv = (double *)v.data;

    Mat rad(rows, cols, CV_64F);
    double *pr = (double *)rad.data;
    Mat k0(rows, cols, CV_32S), k1(rows, cols, CV_32S), f(rows, cols, CV_64F);
    int *pk = (int *)k0.data, *pk1 = (int *)k1.data;
    int kstep = k0.step / sizeof(int), koffset;
    
    double *pf = (double *)f.data, fk;

    Mat wheel = makeColorWheel();
    double *pw = (double *)wheel.data;
    int wstep = wheel.step / sizeof(double), ncols = wheel.rows;

    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c;

            // rad = sqrt(u^2 + v^2)
            pr[offset] = pv[offset]*pv[offset] + pu[offset]*pu[offset];
            pr[offset] = sqrt(pr[offset]);
            
            // a = atan2(-v, -u) / pi
            fk = atan2(-pv[offset], -pu[offset]) / PI;

            // -1~1 mapped to 1~ncols
            fk = (fk+1) / 2 * (ncols-1) + 1;

            // 1, 2, ..., ncols
            koffset = r * kstep + c;
            pk[koffset] = (int)floor(fk);

            if (pk[koffset] == ncols)
                pk1[koffset] = 1;
            else
                pk1[koffset] = pk[koffset] + 1;

            pf[koffset] = fk - pk[koffset];
        }
    }

    Mat img(rows, cols, CV_8UC3);
    int istep = img.step / sizeof(uchar), ioffset;
    uchar *pi = img.data;
    double col0, col1, col;
    bool idx;
    
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c;
            koffset = r * kstep + c;
            ioffset = r * istep + c;
            idx = (pr[offset] <= 1);
        
            for (int k = 0; k < wheel.cols; k++)
            {
                col0 = pw[pk[koffset] * wstep + k] / 255;
                col1 = pw[pk1[koffset] * wstep + k] / 255;
                col = (1 - pf[offset]) * col0 + pf[offset] * col1;

                // increase saturation with radius
                if (idx) col = 1 - pr[offset] * (1 - col);
                // out of range
                else col = col * 0.75;

                pi[ioffset + k] = (uchar)floor(255 * col);
            }
        }
    }

    return img;
}

Mat makeColorWheel()
{
    const int RY = 15;
    const int YG = 6;
    const int GC = 4;
    const int CB = 11;
    const int BM = 13;
    const int MR = 6;

    int rows = RY + YG + GC + CB + BM + MR;
    Mat wheel = Mat::zeros(rows, 3, CV_64F);
    int step = wheel.step / sizeof(double);
    double *pw = (double *)wheel.data;

    // RY
    int start = 0;
    for (int r = start; r < start + RY; r++)
    {
        pw[r * step] = 255;
        pw[r * step + 1] = floor(255 * r / RY);
    }
    start += RY;
    
    // YG
    for (int r = start; r < start + YG; r++)
    {
        pw[r * step] = 255 - floor(255 * (r-start) / YG);
        pw[r * step + 1] = 255;
    }
    start += YG;

    // GC
    for (int r = start; r < start + GC; r++)
    {
        pw[r * step + 1] = 255;
        pw[r * step + 2] = floor(255 * (r-start) / GC);
    }
    start += GC;

    // CB
    for (int r = start; r < start + CB; r++)
    {
        pw[r * step + 1] = 255 - floor(255 * (r-start) / CB);
        pw[r * step + 2] = 255;
    }
    start += CB;

    // BM
    for (int r = start; r < start + BM; r++)
    {
        pw[r * step + 2] = 255;
        pw[r * step] = floor(255 * (r-BM) / BM);
    }
    start += BM;
    
    // MR
    for (int r = start; r < start + MR; r++)
    {
        pw[r * step + 2] = 255 - floor(255 * (r-MR) / MR);
        pw[r * step] = 255;
    }
    
    return wheel;
}
