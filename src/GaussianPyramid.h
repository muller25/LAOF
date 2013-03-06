#ifndef _GaussianPyramid_H
#define _GaussianPyramid_H

#include "Image.h"
#include "ImageProcess.h"

class GaussianPyramid
{
public:
    GaussianPyramid():pPyr(NULL){}
    
    virtual ~GaussianPyramid()
    {
        if (pPyr != NULL)
            delete []pPyr;

        pPyr = NULL;
    }

    void build(const DImage &im, double ratio, int minWidth);
    
    inline DImage& operator[](int i){return pPyr[i];}
    inline int nLevels() const {return levels;}
    
private:
    DImage *pPyr;
    int levels;
};

#endif
