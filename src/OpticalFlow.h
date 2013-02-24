#ifndef _OF_H
#define _OF_H

#include <opencv2/core/core.hpp>
#include <cmath>
using namespace cv;

#include "Maths.h"

#define EPSILON 10e-6

class OpticalFlow
{
public:
    inline double psi(double s, double e=EPSILON);
    inline double psi_d(double s, double e=EPSILON);
    inline double phi(double s, double e=EPSILON);
    inline double phi_d(double s, double e=EPSILON);

    template<class TI, class TF>
    Mat warpImage(Mat& im1, Mat &im2, Mat &u, Mat &v);

    template<class TI, class TF>
    Mat warpImage(Mat& im, Mat &u, Mat &v, bool zero_fill=true);
    
    int compute(Mat &im1, Mat &im2, Mat &u, Mat &v, double a_s);
};

// warping im2 with respect to flow [u v]T
// for point moves out of boundary, we keep the original value in im1
// use bilinear interpolation for warping
template<class TI, class TF>
Mat OpticalFlow::warpImage(Mat &im1, Mat &im2, Mat &u, Mat &v)
{
    assert(Maths::matchAll(im1, im2) && Maths::matchAll(u, v) &&
           Maths::match2D(im1, u));

    int r, c, k, ioffset, foffset;
    double nx, ny;
    int rows = im1.rows;
    int cols = im1.cols;
    int channels = im1.channels();
    int istep = im1.step / sizeof(TI);
    int fstep = u.step / sizeof(TF);
    Mat warp = Mat::zeros(rows, cols, im1.type());
    TI *pim1 = (TI *)im1.data, *pwarp = (TI *)warp.data;
    TF *pu = (TF *)u.data, *pv = (TF *)v.data;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            ioffset = c * istep + c * channels;
            foffset = c * fstep + c;
            nx = c + pu[foffset];
            ny = r + pv[foffset];
                
            if (nx < 0 || nx >= cols || ny < 0 || ny >= rows)
            {
                for (k = 0; k < channels; k++)
                    pwarp[ioffset+k] = pim1[ioffset+k];

                continue;
            }

            Maths::BiInterpolate(im2, nx, ny, pwarp+ioffset);
        }
    }
    
    return warp;
}

// warping im with respect to flow [u v]T
// for point moves out of boundary, we keep the original value in im
// use bilinear interpolation for warping
template<class TI, class TF>
Mat OpticalFlow::warpImage(Mat& im, Mat &u, Mat &v, bool zero_fill)
{
    assert(Maths::matchAll(u, v) && Maths::match2D(im, u));

    int r, c, k, ioffset, foffset;
    double nx, ny;
    int rows = im.rows;
    int cols = im.cols;
    int channels = im.channels();
    Mat warp = Mat::zeros(rows, cols, im.type());
    int istep = im.step / sizeof(TI);
    int fstep = u.step / sizeof(TF);
    TI *pim = (TI *)im.data, *pwarp = (TI *)warp.data;
    TF *pu = (TF *)u.data, *pv = (TF *)v.data;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            ioffset = r * istep + c * channels;
            foffset = r * fstep + c;
            nx = c + pu[foffset];
            ny = r + pv[foffset];
            if (nx < 0 || nx >= cols || ny >= rows || ny < 0)
            {
                // set pixel that moves out of boundary to zero
                if (!zero_fill)
                {
                    for (k = 0; k < channels; k++)
                        pwarp[ioffset+k] = pim[ioffset+k];
                }
                
                continue;
            }

            Maths::BiInterpolate(im, nx, ny, pwarp+ioffset);
        }
    }

    return warp;
}

#endif
