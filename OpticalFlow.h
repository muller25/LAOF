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
    Mat warpI1(Mat& im1, Mat &im2, Mat &u, Mat &v);

    Mat warpI1(Mat& im1, Mat &im2, Mat &flow);
    int compute(Mat &im1, Mat &im2, Mat &u, Mat &v, double a_s);
};

template<class TI, class TF>
Mat OpticalFlow::warpI1(Mat &im1, Mat &im2, Mat &u, Mat &v)
{
    assert(Maths::match(im1, im2, true));
    int r, c, k, offset;
    double nx, ny;
    TI *ptr2, *ptrw;
    TF *ptru, *ptrv;
    int rows = im1.rows;
    int cols = im1.cols;
    int channels = im1.channels();
    Mat warpedI1 = Mat::zeros(rows, cols, im1.type());

    for (r = 0; r < rows; r++)
    {
        ptru = u.ptr<TF>(r);
        ptrv = v.ptr<TF>(r);
        ptr2 = im2.ptr<TI>(r);
        ptrw = warpedI1.ptr<TI>(r);
        
        for (c = 0; c < cols; c++)
        {
            nx = c + ptru[c];
            ny = r + ptrv[c];
            offset = c * channels;
            
            if (nx < 0 || nx >= cols || ny < 0 || ny >= rows)
            {
                for (k = 0; k < channels; k++)
                    ptrw[offset+k] = ptr2[offset+k];

                continue;
            }

            // im1 or im2?
            Maths::bilinear_interpolate(im1, nx, ny, ptrw+offset);
        }
    }
    
    return warpedI1;
}

#endif
