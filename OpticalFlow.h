#ifndef _OF_H
#define _OF_H

#include <opencv2/core/core.hpp>
#include <cmath>
using namespace cv;

#define EPSILON 10e-6

class OpticalFlow
{
    inline double psi(double s, double e=EPSILON);
    inline double psi_d(double s, double e=EPSILON);
    inline double phi(double s, double e=EPSILON);
    inline double phi_d(double s, double e=EPSILON);
    int compute(Mat &im1, Mat &im2, Mat &v, Mat &u, double a_s);
};

#endif
