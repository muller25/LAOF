#ifndef _MATHS_H
#define _MATHS_H

#include "Constants.h"

#include <opencv2/core/core.hpp>
using namespace cv;

class Maths
{
public:
    static Mat dx(const Mat &src);
    static Mat dy(const Mat &src);
    static Mat dxx(const Mat &src);
    static Mat dyy(const Mat &src);
    static Mat dxy(const Mat &src);
    static int sor_solver(const Mat &A, const Mat &b, Mat &x,
                          int nIters=100, double w=1.8);
};
    
#endif
