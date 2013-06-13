#ifndef TYPES_H
#define TYPES_H

#include <cv.h>
using cv::Vec;

typedef float PERCISION;
#define CV_PERCSION CV_32F

typedef int LABEL;
#define CV_LABEL CV_32S

typedef Vec<PERCISION, 3> PIXEL;
typedef Vec<PERCISION, 2> MOTION;

#endif
