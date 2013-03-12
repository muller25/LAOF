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

    void estLapNoise(const DImage &im1, const DImage &im2);

    void im2feature(DImage &feature, const DImage &im);

    void genInImageMask(DImage &mask, const DImage &mask1, const DImage &mask2,
                        const DImage &u, const DImage &v);

    // coarse to fine flow estimation
    void c2fFlow(DImage &u, DImage &v, const DImage &im1, const DImage &im2,
                 double as=0.012, double ratio=0.75, int minWidth=20,
                 int nOutIter=7, int nInIter=1, int nSORIter=30);

    void SORSolver(DImage &u, DImage &v, DImage &warp,
                   const DImage &im1, const DImage &im2, 
                   double as, int nOutIter, int nInIter, int nSORIter);

    void biC2FFlow(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                   const DImage &im1, const DImage &im2,
                   const DImage &mask1, const DImage &mask2,
                   double as=0.012, double ap=0.012, double ratio=0.75, int minWidth=20,
                   int nBiIter=7, int nIRLSIter=1, int nSORIter=30);

    void biIRLS(DImage &du, DImage &dv,
                const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                const DImage &u, const DImage &v,
                const DImage &ur, const DImage &vr,
                double as, double ap, int nIRLSIter, int nSORIter);

    void adIRLS(DImage &du, DImage &dv,
                const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                const DImage &u, const DImage &v,
                const DImage &ur, const DImage &vr,
                double as, double ap, int nIRLSIter, int nSORIter);

    void stFlow(std::vector<DImage> &u, std::vector<DImage> &v,
                std::vector<DImage> &ur, std::vector<DImage> &vr,
                std::vector<DImage> &mask, const std::vector<DImage> &im,
                int idx, double as, double ap, int nBiIter, int nIRLSIter, int nSORIter);
    
    void stFlow(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                const DImage &im1, const DImage &im2,
                const DImage &mask1, const DImage &mask2,
                const DImage &pu, const DImage &pv,
                const DImage &nu, const DImage &nv,
                const DImage &pur, const DImage &pvr,
                const DImage &nur, const DImage &nvr,
                double as, double ap, int nBiIter, int nIRLSIter, int nSORIter);
    
    void adIRLS3(DImage &du, DImage &dv,
                 const DImage &Ix, const DImage &Iy, const DImage &It,
                 const DImage &mask, const DImage &pphid,
                 const DImage &pu, const DImage &pv,
                 const DImage &u, const DImage &v,
                 const DImage &nu, const DImage &nv,
                 const DImage &ur, const DImage &vr,
                 double as, double ap, int nIRLSIter, int nSORIter);

private:
    std::vector<double> lapPara;
};

#endif
