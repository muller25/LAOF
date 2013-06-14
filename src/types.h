#ifndef TYPES_H
#define TYPES_H

#include <cv.h>
using cv::Vec;

typedef uchar IMAGE;
#define CV_IMAGE CV_8U

typedef float PERCISION;
#define CV_PERCISION CV_32F

typedef int LABEL;
#define CV_LABEL CV_32S

typedef Vec<IMAGE, 3> PIXEL;
typedef Vec<PERCISION, 2> MOTION;

#endif
