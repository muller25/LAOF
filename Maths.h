#ifndef _MATHS_H
#define _MATHS_H

#include "Constants.h"

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
#include <cmath>
using namespace std;

class Maths
{
public:
    static Mat dx(const Mat &src);
    static Mat dy(const Mat &src);
    static Mat dxx(const Mat &src);
    static Mat dyy(const Mat &src);
    static Mat dxy(const Mat &src);
    static int sor_solver(const Mat &A, const Mat &b, Mat &x,
                          int nIters=100, double w=1.8);
    static bool match(const Mat& m1, const Mat &m2, bool depth=false);

    template<class TF, class TW>
    static Mat weighted_laplacian(const Mat &flow, const Mat &weight);

    template<class TF, class TW>
    static Mat weighted_laplacian3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
                                    const Mat &weight, const Mat &nweight);
    
    template<class T>
    static void bilinear_interpolate(Mat &im, double x, double y, T *res);
};

/*
  [description]
  first apply backward difference, then apply forward difference
  to approximate weighted laplacian
  [params]
  flow - flow matrix, it may contain multi-channel (in)
  weight - weight matrix, it only contains one-channel (in)
  [return]
  weighted laplacian matrix
 */
template<class TF, class TW>
Mat Maths::weighted_laplacian(const Mat &flow, const Mat &weight)
{
    int r, c, k, offset, woffset, loffset;
    int rows = flow.rows;
    int cols = flow.cols;
    int channels = flow.channels();
    int step = flow.step / sizeof(TF);
    TF *fptr = (TF *)flow.data;

    int wstep = weight.step / sizeof(TW);
    TW *wptr = (TW *)weight.data;
    
    Mat lap = Mat::zeros(rows, cols, CV_64FC(channels));
    int lstep = lap.step / sizeof(double);
    double *lptr = (double *)lap.data;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            woffset = r * wstep + c;
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                loffset = r * lstep + c * channels + k;
                if (c < cols-1)
                    lptr[loffset] += wptr[woffset+1] * (fptr[offset+channels] - fptr[offset]);
                if (c > 0)
                    lptr[loffset] -= wptr[woffset] * (fptr[offset] - fptr[offset-channels]);
                if (r < rows-1)
                    lptr[loffset] += wptr[woffset+wstep] * (fptr[offset+step] - fptr[offset]);
                if (r > 0)
                    lptr[loffset] -= wptr[woffset] * (fptr[offset] - fptr[offset-step]);
            }
        }
    }

    return lap;
}

/*
 */
template<class TF, class TW>
Mat Maths::weighted_laplacian3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
                                const Mat &weight, const Mat &nweight)
{
    // assert(match(pflow, flow, true) && match(flow, nflow, true) && flow.channels() == 2 &&
    //        match(weight, nweight, true) && weight.channels() == 1 && 
    //        flow.size() == weight.size() &&
    //        weight.depth() == CV_64F && flow.depth() == CV_64F);
    
    int rows = flow.rows;
    int cols = flow.cols;
    int channels = flow.channels();
    int step = 0;//flow.step / get_step(flow.depth());
    int wstep = 0;//weight.step / get_step(weight.depth());
    int r, c, k, offset, nr, nc, noffset, woffset;
    double *pfptr, *fptr, *nfptr, *wptr, *nwptr, *l3ptr;

    Mat lap3d = weighted_laplacian<double, double>(flow, weight);
    l3ptr = (double *)lap3d.data;
    pfptr = (double *)pflow.data;
    fptr = (double *)flow.data;
    nfptr = (double *)nflow.data;
    wptr = (double *)weight.data;
    nwptr = (double *)nweight.data;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            nr = r + pfptr[offset];
            nc = c + pfptr[offset+1];
            if (nr < 0 || nr >= rows || nc >= cols || nc < 0) continue;

            woffset = r * wstep + c;
            for (k = 0; k < channels; k++)
            {
                noffset = nr * step + nc * channels + k;
                l3ptr[noffset] -= wptr[woffset] * (fptr[noffset] - pfptr[offset+k]);
            }
        }
    }

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            nr = fptr[offset] + r;
            nc = fptr[offset+1] + c;
            if (nr < 0 || nr >= rows || nc >= cols || nc < 0) continue;

            noffset = nr * step + nc * channels;
            woffset = nr * wstep + nc;
            for (k = 0; k < channels; k++)
                l3ptr[offset+k] += nwptr[woffset] * (nfptr[noffset+k] - fptr[offset+k]);
        }
    }
    
    return lap3d;
}

template<class T>
void Maths::bilinear_interpolate(Mat &im, double x, double y, T *res)
{
    int r, c, k, nx, ny;
    T *ptr;
    int ix = x, iy = y;
    double dx = x-ix, dy = y-iy, s;
    int channels = im.channels();
    int rows = im.rows;
    int cols = im.cols;
    
    for (r = 0; r <= 1; r++)
    {
        ny = min(max(iy+r, 0), rows-1);
        ptr = im.ptr<T>(ny);
        for (c = 0; c <= 1; c++)
        {
            nx = min(max(ix+c, 0), cols-1);
            s = fabs(1-r-dy) * fabs(1-c-dx);
            for (k = 0; k < channels; k++)
                res[k] += ptr[nx*channels+k] * s;
        }
    }
}

#endif
