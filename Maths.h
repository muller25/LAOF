#ifndef _MATHS_H
#define _MATHS_H

#include "Constants.h"

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
using namespace std;

class Maths
{
public:
    static Mat dx(const Mat &src);
    static Mat dy(const Mat &src);
    static Mat dxx(const Mat &src);
    static Mat dyy(const Mat &src);
    static Mat dxy(const Mat &src);
    static Mat laplace3D(const Mat &im0, const Mat &im1, const Mat &im2);
    static int sor_solver(const Mat &A, const Mat &b, Mat &x,
                          int nIters=100, double w=1.8);
};
    
#endif
