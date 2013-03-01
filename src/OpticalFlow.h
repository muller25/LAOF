#ifndef _OF_H
#define _OF_H

#include "Maths.h"
#include "GaussianPyramid.h"

#include <cv.h>
using namespace cv;

#include <cmath>
#include <cstdio>
#include <vector>

#define EPSILON 1e-6

class OpticalFlow
{
public:
    OpticalFlow()
    {
        lapPara.resize(9, 0.02);
    }
    virtual ~OpticalFlow()
    {
        lapPara.clear();
    }
    
    double psi(const double s, const double e=EPSILON)
    {
        return sqrt(s + e);
    }
    
    double psi_d(const double s, const double e=EPSILON)
    {
        return 0.5 / sqrt(s + e);
    }

    void psi_d(const Mat &Ix, const Mat &Iy, const Mat &Iz,
               const Mat &du, const Mat &dv, Mat &res);
    
    double phi(const double s, const double e=EPSILON)
    {
        return sqrt(s + e);
    }
    
    double phi_d(const double s, const double e=EPSILON)
    {
        return 0.5 / sqrt(s + e);
    }

    void phi_d(const Mat &u, const Mat &v, Mat &res);
    
    void warpImage(const Mat& im1, const Mat &im2,
                   const Mat &u, const Mat &v, Mat &warp);
    
    void getGrads(const Mat &im1, const Mat &im2, Mat &Ix, Mat &Iy, Mat &Iz);
    
    void compute(const Mat &im1, const Mat &im2, Mat &warp,
                 Mat &u, Mat &v, const double a_s,
                 const int nOutIter=100, const int nInIter=100, const int nSORIter=100);

    void sanityCheck(const Mat &Ix, const Mat &Iy, const Mat &It,
                     const Mat &du, const Mat &dv);

    void estLapNoise(const Mat &im1, const Mat &im2);

    void Coarse2FineFlow(const Mat &im1, const Mat &im2, Mat &warp, Mat &u, Mat &v,
                         const double a_s, const double ratio, const int minWidth,
                         const int nOutIter, const int nInIter, const int nSORIter);

    void im2feature(const Mat &im, Mat &feature);

    template<class T>
    void threshold(Mat &im);
    
private:
    std::vector<double> lapPara;
};

template<class T>
void OpticalFlow::threshold(Mat &im)
{
    int rows = im.rows, cols = im.cols, channels = im.channels();
    int step = im.step / sizeof(double), offset;
    double *p = (double *)im.data;
    const T rangeMin = 0;
    T rangeMax;
    if (typeid(T) == typeid(float) || typeid(T) == typeid(double) || typeid(T) == typeid(long double))
        rangeMax = 1;
    else
        rangeMax = 255;
    
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            for (int k = 0; k < channels; k++)
                p[offset+k] = std::min(std::max(p[offset+k], rangeMin), rangeMax);
        }
    }
    
}
#endif
