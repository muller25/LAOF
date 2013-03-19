#ifndef _OF_H
#define _OF_H

#include "Image.h"

class OpticalFlow
{
public:
    inline void psi_d(DImage &res,
                      const DImage &Ix, const DImage &Iy, const DImage &It,
                      const DImage &du, const DImage &dv);

    inline void phi_d(DImage &res, const DImage &u, const DImage &v);
    
    void getGrads(DImage &Ix, DImage &Iy, DImage &It,
                  const DImage &im1, const DImage &im2);

    void im2feature(DImage &feature, const DImage &im);

    void genInImageMask(DImage &mask, const DImage &mask1, const DImage &mask2,
                        const DImage &u, const DImage &v);

    // coarse to fine strategy for normal optical flow
    void c2fFlow(DImage &u, DImage &v, const DImage &im1, const DImage &im2,
                 double as=0.012, double ratio=0.75, int minWidth=20,
                 int nOutIter=7, int nInIter=1, int nSORIter=30);

    void SORSolver(DImage &u, DImage &v, DImage &warp,
                   const DImage &im1, const DImage &im2, 
                   double as, int nOutIter, int nInIter, int nSORIter);

    // coarse to fine strategy for biDir optical flow
    void biC2FFlow(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                   const DImage &im1, const DImage &im2,
                   const DImage &mask1, const DImage &mask2,
                   double as=0.012, double ap=0.012, double ratio=0.75, int minWidth=20,
                   int nBiIter=7, int nIRLSIter=1, int nSORIter=30);

    // biDir optical flow
    void biIRLS(DImage &du, DImage &dv,
                const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                const DImage &u, const DImage &v,
                const DImage &ur, const DImage &vr,
                double as, double ap, int nIRLSIter, int nSORIter);

    // adaptive optical flow based on biDir flow
    void adIRLS(DImage &du, DImage &dv,
                const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                const DImage &u, const DImage &v,
                const DImage &ur, const DImage &vr,
                double as, double ap, int nIRLSIter, int nSORIter);

    // coarse to fine strategy for adaptive flow based on motion covariance
    void adC2FFlow(DImage &u, DImage &v, 
                   const DImage &im1, const DImage &im2,
                   const DImage &mask1, const DImage &mask2,
                   double as, double ratio, int minWidth,
                   int nOutIter, int nIRLSIter, int nSORIter);

    // adpative optical flow based on motion covariance
    void adIRLS2(DImage &du, DImage &dv, const DImage &D,
                 const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                 const DImage &u, const DImage &v,
                 double as, int nIRLSIter, int nSORIter);

    // coarse to fine strategy for adaptive spatial temporal smooth optical flow
    // based on motion covariance
    void stC2FFlow(std::vector<DImage> &u, std::vector<DImage> &v, 
                   const std::vector<DImage> &im,
                   const std::vector<DImage> &mask, int i0,
                   double as, double at, double ratio, int minWidth,
                   int nOutIter, int nIRLSIter, int nSORIter);

    // temporal smooth operator for optical flow
    // based on motion covariance
    void temporalSmooth(DImage &u0, DImage &v0,
                        const DImage &u1, const DImage &v1,
                        const DImage &im0, const DImage &im1,
                        const DImage &mask0, const DImage &mask1,
                        double as, double at,
                        int nOutIter, int nIRLSIter, int nSORIter);
};

#endif
