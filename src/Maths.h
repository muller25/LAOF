#ifndef _MATHS_H
#define _MATHS_H

#include <cv.h>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <cmath>
using namespace std;

#define ESP 10e-6

Mat dx(const Mat &src);
Mat dy(const Mat &src);
Mat dxx(const Mat &src);
Mat dyy(const Mat &src);
Mat dxy(const Mat &src);
    
bool match2D(const Mat& m1, const Mat &m2);
bool match3D(const Mat& m1, const Mat &m2);
bool matchAll(const Mat& m1, const Mat &m2);

void weighted_lap(const Mat &flow, const Mat &weight, Mat &dst);

void collapse(Mat &src, Mat &dst);

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
void biInterpolate(Mat &im, double x, double y, T *res)
{
    int r, c, k, nx, ny, offset;
    int ix = x, iy = y;
    double s, dx = x-ix, dy = y-iy;
    int rows = im.rows, cols = im.cols, channels = im.channels();
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
