#ifndef _OF_H
#define _OF_H

#include <opencv2/core/core.hpp>
using namespace cv;

class OpticalFlow
{
public:
    static int SOR_solver(const Mat &A, const Mat &b, Mat &x,
                          int nIters, double w = 1.8);
};

#endif
