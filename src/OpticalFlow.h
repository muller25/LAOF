#ifndef _OF_H
#define _OF_H

#include "Maths.h"

#include <cv.h>
using namespace cv;

#include <cmath>
#include <iostream>
using namespace std;

#define EPSILON 10e-6

class OpticalFlow
{
public:
    double psi(double s, double e=EPSILON){return sqrt(s + e);}
    double psi_d(double s, double e=EPSILON){return 0.5 / sqrt(s + e);}
    void psi_d(Mat &Ix, Mat &Iy, Mat &Iz, Mat &du, Mat &dv, Mat &res);
    
    double phi(double s, double e=EPSILON){return sqrt(s + e);}
    double phi_d(double s, double e=EPSILON){return 0.5 / sqrt(s + e);}
    void phi_d(Mat &u, Mat &v, Mat &res);
    
    void warpImage(Mat& im1, Mat &im2, Mat &u, Mat &v, Mat &warp);
    void getGrads(Mat &im1, Mat &im2, Mat &Ix, Mat &Iy, Mat &Iz);
    
    void compute(Mat &im1, Mat &im2, Mat &warp, Mat &u, Mat &v, double a_s,
                 int nInIter=100, int nOutIter=100, int nSORIter=1000);
};

#endif
