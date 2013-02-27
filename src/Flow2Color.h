#ifndef _FLOW2COLOR_H
#define _FLOW2COLOR_H

#include <cstdio>

#include <cv.h>
using namespace cv;

#include "Maths.h"

void flow2color(Mat &u, Mat &v, Mat &flow, Mat &idx);
Mat computeColor(Mat &u, Mat &v);
Mat makeColorWheel();

#endif
