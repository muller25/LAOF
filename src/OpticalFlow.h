#ifndef _OF_H
#define _OF_H

#include "Maths.h"

#include <cv.h>
using namespace cv;

#include <cmath>
#include <cstdio>

#define EPSILON 1e-6

class OpticalFlow
{
public:
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

};

#endif
