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
    static int SORSolver(const Mat &A, const Mat &b, Mat &x,
                         int nIters=100, double w=1.8);
    
    static bool match2D(const Mat& m1, const Mat &m2);
    static bool match3D(const Mat& m1, const Mat &m2);
    static bool matchAll(const Mat& m1, const Mat &m2);

    template<class TF, class TW>
    static Mat weighted_lap(const Mat &flow, const Mat &weight);

    template<class TF, class TW>
    static Mat weighted_lap3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
                              const Mat &weight, const Mat &nweight);
    
    template<class T>
    static void BiInterpolate(Mat &im, double x, double y, T *res);
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
Mat Maths::weighted_lap(const Mat &flow, const Mat &weight)
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
Mat Maths::weighted_lap3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
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

    Mat lap3d = weighted_lap<double, double>(flow, weight);
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


/*
  [description]
  计算二维矩阵im中位置(x,y)的双线性插值结果,并保存到res
  [params]
  im - 二维矩阵，可以多通道
  x - 目标位置的横坐标
  y - 目标位置的竖坐标
  res - 插值结果，通道数应该和im一致
  [return]
  无
*/
template<class T>
void Maths::BiInterpolate(Mat &im, double x, double y, T *res)
{
    assert(res != NULL);
    
    int r, c, k, nx, ny, offset;
    int ix = x, iy = y;
    double s, dx = x-ix, dy = y-iy;
    int channels = im.channels();
    int rows = im.rows;
    int cols = im.cols;
    int step = im.step / sizeof(T);
    T *ptr = (T *)im.data;
    double tmp[3] = {0};
    
    for (r = 0; r <= 1; r++)
    {
        ny = min(max(iy+r, 0), rows-1);
        for (c = 0; c <= 1; c++)
        {
            nx = min(max(ix+c, 0), cols-1);
            s = fabs(1-r-dy) * fabs(1-c-dx);
            offset = nx * channels + ny * step;

            for (k = 0; k < channels; k++)
                tmp[k] += ptr[offset+k] * s;
        }
    }

    for (k = 0; k < channels; k++)
        res[k] += tmp[k];
}

#endif
