#ifndef _GaussianPyramid_H
#define _GaussianPyramid_H

#include <cv.h>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

class GaussianPyramid
{
public:
    GaussianPyramid()
    {
        pPyr = NULL;
    }

    virtual ~GaussianPyramid()
    {
        if (pPyr != NULL)
            delete []pPyr;

        pPyr = NULL;
    }
    
    void ConstructPyramid(const Mat &im, double ratio, int minWidth);
    void ConstructPyramid(const Mat &im, int minWidth);
    
    Mat& operator[](int i)
    {
        return pPyr[i];
    }

    int nLevels()
    {
        return levels;
    }
    
private:
    Mat *pPyr;
    int levels;
};

#endif
