#ifndef _OF_H
#define _OF_H

#include "Image.h"
#include "Maths.h"
#include "ImageProcess.h"
#include "GaussianPyramid.h"

#include <cmath>
#include <cstdio>
#include <vector>

class OpticalFlow
{
public:
    OpticalFlow(){lapPara.resize(9, 0.02);}
    virtual ~OpticalFlow(){lapPara.clear();}

    inline void psi_d(DImage &res,
                      const DImage &Ix, const DImage &Iy, const DImage &It,
                      const DImage &du, const DImage &dv);

    inline void phi_d(DImage &res, const DImage &u, const DImage &v);
    
    void warpImage(DImage &warp, const DImage &im1, const DImage &im2,
                   const DImage &u, const DImage &v);
    
    void getGrads(DImage &Ix, DImage &Iy, DImage &It,
                  const DImage &im1, const DImage &im2);
    
    void SORSolver(DImage &u, DImage &v, DImage &warp,
                   const DImage &im1, const DImage &im2, 
                   double a_s, int nOutIter=7, int nInIter=1, int nSORIter=30);

    void estLapNoise(const DImage &im1, const DImage &im2);

    // coarse to fine flow estimation
    void c2fFlow(DImage &u, DImage &v, const DImage &im1, const DImage &im2,
                 double a_s, double ratio, int minWidth,
                 int nOutIter, int nInIter, int nSORIter);

    void im2feature(DImage &feature, const DImage &im);

private:
    std::vector<double> lapPara;
};

#endif
