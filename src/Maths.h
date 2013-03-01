#ifndef _MATHS_H
#define _MATHS_H

#include <cv.h>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <cmath>
using namespace std;

#define ESP 1e-6

Mat dx(const Mat &src);
Mat dy(const Mat &src);
Mat dxx(const Mat &src);
Mat dyy(const Mat &src);
Mat dxy(const Mat &src);
    
bool match2D(const Mat& m1, const Mat &m2);
bool match3D(const Mat& m1, const Mat &m2);
bool matchAll(const Mat& m1, const Mat &m2);

void weighted_lap(const Mat &flow, const Mat &weight, Mat &dst);

template<class T>
void collapse(const Mat &src, Mat &dst)
{
    int rows = src.rows, cols = src.cols, channels = src.channels();
    T *ps = (T *)src.data;
    int sstep = src.step / sizeof(T), soffset;

    if (channels == 1)
    {
        src.copyTo(dst);
        return;
    }

    dst.create(rows, cols, src.depth());
    dst.setTo(0);
    T  *pd = (T *)dst.data;
    int dstep = dst.step / sizeof(T), doffset, r, c, k;
    double tmp;
        
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            soffset = r * sstep + c * channels;
            doffset = r * dstep + c;
            for (tmp = 0, k = 0; k < channels; k++)
                tmp += ps[soffset + k];
            
            pd[doffset] = (T)(tmp / channels);
        }
    }
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
inline void biInterpolate(const Mat &im, const double x, const double y, T *res)
{
    int r, c, k, nx, ny, offset;
    int ix = x, iy = y;
    double s, dx, dy;
    int rows = im.rows, cols = im.cols, channels = im.channels();
    int step = im.step / sizeof(T);
    T *ptr = (T *)im.data;
    double tmp[3] = {0};

    dx = std::max(std::min(x-floor(x), 1.), 0.);
    dy = std::max(std::min(y-floor(y), 1.), 0.);
    for (r = 0; r <= 1; r++)
    {
        ny = std::min(std::max(iy+r, 0), rows-1);
        for (c = 0; c <= 1; c++)
        {
            nx = std::min(std::max(ix+c, 0), cols-1);
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
