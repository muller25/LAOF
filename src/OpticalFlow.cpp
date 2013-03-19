#include "OpticalFlow.h"
#include "Maths.h"
#include "ImageProcess.h"
#include "GaussianPyramid.h"

#include "ImageIO.h"
#include "Flow2Color.h"

#include <cmath>
#include <cstdio>
#include <vector>

// calculate res = psi_1st((It + Ix * du + Iy * dv)^2)
void OpticalFlow::psi_d(DImage &res, const DImage &Ix, const DImage &Iy,
                        const DImage &It, const DImage &du, const DImage &dv)
{
    assert(Ix.match3D(Iy) && Iy.match3D(It) && du.match3D(dv) &&
           It.match2D(dv) && du.nChannels() == 1);

    const double epsilon = 1e-6;
    
    int width = Ix.nWidth(), height = Ix.nHeight(), channels = Ix.nChannels();
    int io, fo;
    double tmp;
    
    res.create(width, height, channels);    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            fo = h * width + w;
            io = fo * channels;
            for (int k = 0; k < channels; ++k)
            {
                 tmp = It[io+k] + Ix[io+k]*du[fo] + Iy[io+k]*dv[fo];
                tmp *= tmp;
                res[io+k] = 0.5 / sqrt(tmp + epsilon);
            }
        }
    }
}

// calculate res = phi_1st(u^2 + v^2)
void OpticalFlow::phi_d(DImage &res, const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const double epsilon = 1e-6;
    DImage ux, uy, vx, vy;
    double tmp;    

    grad1st(ux, uy, u);
    grad1st(vx, vy, v);
    res.create(u.nWidth(), u.nHeight());
    
    for (int i = 0; i < res.nElements(); i++)
    {
        tmp = ux[i]*ux[i] + uy[i]*uy[i] + vx[i]*vx[i] + vy[i]*vy[i];
        res[i] = 0.5 / sqrt(tmp + epsilon);
    }
}

// get Ix, Iy, It
// smooth im1 -> sim1, im2 -> sim2
// im = 0.4 * sim1 + 0.6 * sim2
// Ix = dx(im), Iy = dy(im), It = sim2 - sim1
void OpticalFlow::getGrads(DImage &Ix, DImage &Iy, DImage &It,
                           const DImage &im1, const DImage &im2)
{
    // double filter[5] = {0.02, 0.11, 0.74, 0.11, 0.02};
    // DImage sIm1, sIm2, im;

    // filtering(sIm1, im1, filter, 2, filter, 2);
    // filtering(sIm2, im2, filter, 2, filter, 2);
    // addWeighted(im, sIm1, 0.4, sIm2, 0.6);
    // grad1st(Ix, Iy, im);
    // substract(It, sIm2, sIm1);

    // suppose we use Guided Filter before, hence we don't use any filter here
    DImage im;
    addWeighted(im, im1, 0.4, im2, 0.6);
    grad1st(Ix, Iy, im);
    substract(It, im2, im1);
}

void OpticalFlow::im2feature(DImage &feature, const DImage &im)
{
    int height = im.nHeight(), width = im.nWidth(), channels = im.nChannels();
    const int nchannels = channels + 2;
    const double gamma = 12./9;
    DImage gx, gy;    
    
    feature.create(width, height, nchannels);
    double *pf = feature.ptr();
    int offset, foffset;
    
    if (channels == 1)
    {
        grad1st(gx, gy, im);
        
        // mix channels
        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; w++)
            {
                offset = h * width + w;
                foffset = offset * nchannels;
                pf[foffset] = im[offset];
                pf[foffset + 1] = gx[offset] * gamma;
                pf[foffset + 2] = gy[offset] * gamma;
            }
        }
        
        return;
    }

    if (channels == 3)
    {
        DImage gray;

        desuarate(gray, im);
        grad1st(gx, gy, gray);
        
        // mix channels
        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; ++w)
            {
                offset = h * width + w;
                foffset = offset * nchannels;
                pf[foffset] = gray[offset];
                pf[foffset + 1] = gx[offset] * gamma;
                pf[foffset + 2] = gy[offset] * gamma;
                
                // green - blue
                pf[foffset + 3] = im[offset*3+1] - im[offset*3];
                // green - red
                pf[foffset + 4] = im[offset*3+1] - im[offset*3+2];
            }
        }
        
        return;
    }
    
    im.copyTo(feature);
}

// generate mask for pixels in mask1, moves in image boundary and mask2
void OpticalFlow::genInImageMask(DImage &mask, const DImage &mask1, const DImage &mask2,
                                 const DImage &u, const DImage &v)
{
    assert(mask1.match3D(mask2) && u.match3D(v) && u.match2D(mask1));
    int width = mask1.nWidth(), height = mask1.nHeight();
    double nu, nv;
    int x, y, offset;
    
    mask.create(width, height);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;

            // mask1[offset] == 0, pixel not in mask1
            if (fabs(mask1[offset]) < ESP) continue;
            
            nu = w + u[offset];
            nv = h + v[offset];
            if (nu < 0 || nu > width-1 || nv < 0 || nv > height-1)
            {
                // give another chance to pixel moves out of boundary
                mask[offset] = 1;
                continue;
            }

            // pixel in mask2
            x = nu + 0.5, y = nv + 0.5;
            if (x < 0 || x > width-1 || y < 0 || y > height-1)
                continue;
            
            if (fabs(mask2[y*width+x]) > ESP)
                mask[offset] = mask1[offset] * mask2[y*width+x];
        }
    }
}

// coarse to fine strategy for normal optical flow
void OpticalFlow::c2fFlow(DImage &u, DImage &v, const DImage &im1, const DImage &im2,
                          double as, double ratio, int minWidth,
                          int nOutIter, int nInIter, int nSORIter)
{
    GaussianPyramid pyr1, pyr2;

    pyr1.build(im1, ratio, minWidth);
    pyr2.build(im2, ratio, minWidth);

    DImage tmp, fIm1, fIm2, warp;
    int width, height;
    double xfactor, yfactor;
    
    // iterate from the top level to the bottom
    for (int k = pyr1.nLevels()-1; k >= 0; --k)
    {
        printf("Pyramid level %d\n", k);

        width = pyr1[k].nWidth();
        height = pyr1[k].nHeight();
        im2feature(fIm1, pyr1[k]);
        im2feature(fIm2, pyr2[k]);
        
        // if at the top level
        if (k == pyr1.nLevels()-1)
        {
            u.create(width, height);
            v.create(width, height);
            fIm2.copyTo(warp);
        } else {
            xfactor = width / pyr1[k+1].nWidth();
            yfactor = height / pyr1[k+1].nHeight();

            imresize(tmp, u, width, height);
            multiply(u, tmp, xfactor);

            imresize(tmp, v, width, height);
            multiply(v, tmp, yfactor);

            warpImage(warp, fIm1, fIm2, u, v);
        }

        SORSolver(u, v, warp, fIm1, fIm2, as, nOutIter+k, nInIter, nSORIter+k*3);
        printf("u: %.6f .. %.6f, v: %.6f .. %.6f\n", u.min(), u.max(), v.min(), v.max());
    }
}

// 计算im1到im2的光流，u, v为初始化值，并将结果保存在u,v。a_s是smoothness term的权重
// u, v, warp需要预先初始化
// normal optical flow
void OpticalFlow::SORSolver(DImage &u, DImage &v, DImage &warp,
                            const DImage &im1, const DImage &im2,
                            double as, int nOutIter, int nInIter, int nSORIter)
{
    assert(im1.match3D(im2) && u.match3D(v) && warp.match3D(im1));

    int width = im1.nWidth(), height = im1.nHeight();
    DImage du(width, height), dv(width, height);
    DImage It, Ix, Iy, Ixy, Ix2, Iy2, Ixt, Iyt, uu, vv;
    DImage psi_1st, phi_1st, lapU, lapV;
    DImage ixy, ix2, iy2, iyt, ixt;
    
    // outer fixed point iteration for u, v
    for (int oiter = 0; oiter < nOutIter; oiter++)
    {
        // compute Ix, Iy and It
        getGrads(Ix, Iy, It, im1, warp);

        // inner fixed point iteration for du, dv
        du.set(0);
        dv.set(0);
        for (int iiter = 0; iiter < nInIter; iiter++)
        {
            add(uu, u, du);// uu = u + du
            add(vv, v, dv);// vv = v + dv
            phi_d(phi_1st, uu, vv);
            psi_d(psi_1st, Ix, Iy, It, du, dv);
            
            // components for linear system
            multiply(Ix2, Ix, Ix, psi_1st);
            collapse(ix2, Ix2);

            multiply(Iy2, Iy, Iy, psi_1st);
            collapse(iy2, Iy2);
            
            multiply(Ixy, Ix, Iy, psi_1st);
            collapse(ixy, Ixy);
            
            multiply(Ixt, It, Ix, psi_1st);
            collapse(ixt, Ixt);
            
            multiply(Iyt, It, Iy, psi_1st);
            collapse(iyt, Iyt);
            
            weighted_lap(lapU, u, phi_1st);
            weighted_lap(lapV, v, phi_1st);
            for (int i = 0; i < iyt.nElements(); ++i)
            {
                ixt[i] = -ixt[i] + as * lapU[i];
                iyt[i] = -iyt[i] + as * lapV[i];
            }
            
            // SOR iteration
            du.set(0);
            dv.set(0);
            const double omega = 1.8;
            double l, l_du, l_dv;
            int offset, tmp;

            for (int siter = 0; siter < nSORIter; siter++)
            {
                for (int h = 0; h < height; ++h)
                {
                    for (int w = 0; w < width; ++w)
                    {
                        offset = h * width + w;
                        l_du = 0, l_dv = 0, l = 0;
                        
                        if (h > 0)
                        {
                            tmp = offset - width;
                            l_du += phi_1st[tmp] * du[tmp];
                            l_dv += phi_1st[tmp] * dv[tmp];
                            l -= phi_1st[tmp];
                        }
                        if (h < height-1)
                        {
                            tmp = offset + width;
                            l_du += phi_1st[offset] * du[tmp];
                            l_dv += phi_1st[offset] * dv[tmp];
                            l -= phi_1st[offset];
                        }
                        if (w > 0)
                        {
                            tmp = offset - 1;
                            l_du += phi_1st[tmp] * du[tmp];
                            l_dv += phi_1st[tmp] * dv[tmp];
                            l -= phi_1st[tmp];
                        }
                        if (w < width-1)
                        {
                            tmp = offset + 1;
                            l_du += phi_1st[offset] * du[tmp];
                            l_dv += phi_1st[offset] * dv[tmp];
                            l -= phi_1st[offset];
                        }

                        l *= as;
                        l_du *= as;
                        l_dv *= as;
                        
                        // du
                        l_du = ixt[offset] + l_du - ixy[offset] * dv[offset];
                        du[offset] = (1-omega)*du[offset]+omega/(ix2[offset]-l)*l_du;
                        
                        // dv
                        l_dv = iyt[offset] + l_dv - ixy[offset] * du[offset];
                        dv[offset] = (1-omega)*dv[offset]+omega/(iy2[offset]-l)*l_dv;
                    }
                }

            }
//            printf("du: %.6f .. %.6f, dv: %.6f .. %.6f\n", du.min(), du.max(), dv.min(), dv.max());
        }

        add(u, du);// u += du
        add(v, dv);// v += dv
        warpImage(warp, im1, im2, u, v);
    }
}

// coarse to fine strategy for biDir optical flow
void OpticalFlow::biC2FFlow(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                            const DImage &im1, const DImage &im2,
                            const DImage &mask1, const DImage &mask2,
                            double as, double ap, double ratio, int minWidth,
                            int nBiIter, int nIRLSIter, int nSORIter)
{
    GaussianPyramid pyr1, pyr2, mpyr1, mpyr2;

    pyr1.build(im1, ratio, minWidth);
    pyr2.build(im2, ratio, minWidth);

    // do not smooth while build pyramid
    mpyr1.build(mask1, ratio, minWidth, false);
    mpyr2.build(mask2, ratio, minWidth, false);
    // printf("pyramid constructed\n");
    
    DImage fIm1, fIm2, warpI1, warpI2, Ix, Iy, It, mask;
    DImage du1, dv1, du2, dv2, tmp;
    int width, height;
    double xfactor, yfactor;

    // iterate from the top level to the bottom
    for (int k = pyr1.nLevels()-1; k >= 0; --k)
    {
        // printf("Pyramid level %d\n", k);

        width = pyr1[k].nWidth();
        height = pyr1[k].nHeight();
        im2feature(fIm1, pyr1[k]);
        im2feature(fIm2, pyr2[k]);
        
        // if at the top level
        if (k == pyr1.nLevels()-1)
        {
            u1.create(width, height);
            v1.create(width, height);
            fIm2.copyTo(warpI2);
            
            u2.create(width, height);
            v2.create(width, height);
            fIm1.copyTo(warpI1);
        } else {
            xfactor = (double)width / pyr1[k+1].nWidth();
            yfactor = (double)height / pyr1[k+1].nHeight();
            
            imresize(tmp, u1, width, height);
            multiply(u1, tmp, xfactor);
            imresize(tmp, v1, width, height);
            multiply(v1, tmp, yfactor);
            warpImage(warpI2, fIm1, fIm2, u1, v1);
            
            imresize(tmp, u2, width, height);
            multiply(u2, tmp, xfactor);
            imresize(tmp, v2, width, height);
            multiply(v2, tmp, yfactor);
            warpImage(warpI1, fIm2, fIm1, u2, v2);
        }
        
        mask.create(width, height, 1, 1);        
        for (int l = 0; l < nBiIter; l++)
        {
            // forward flow
            getGrads(Ix, Iy, It, fIm1, warpI2);
            if (l >= nBiIter-2)
                genInImageMask(mask, mpyr1[k], mpyr2[k], u1, v1);
            
//            biIRLS(du1, dv1, Ix, Iy, It, mask, u1, v1, u2, v2, as, ap, nIRLSIter, nSORIter);
            adIRLS(du1, dv1, Ix, Iy, It, mask, u1, v1, u2, v2, as, ap, nIRLSIter, nSORIter);

            // backward flow
            getGrads(Ix, Iy, It, fIm2, warpI1);
            if (l >= nBiIter-2)
                genInImageMask(mask, mpyr2[k], mpyr1[k], u2, v2);
            
//            biIRLS(du2, dv2, Ix, Iy, It, mask, u2, v2, u1, v1, as*pow(ratio, k), ap, nIRLSIter, nSORIter);
            adIRLS(du2, dv2, Ix, Iy, It, mask, u2, v2, u1, v1, as, ap, nIRLSIter, nSORIter);

            add(u1, du1);
            add(v1, dv1);
            add(u2, du2);
            add(v2, dv2);

            warpImage(warpI2, fIm1, fIm2, u1, v1);
            warpImage(warpI1, fIm2, fIm1, u2, v2);

            // printf("u1: %.6f .. %.6f, v1: %.6f .. %.6f\n", u1.min(), u1.max(), v1.min(), v1.max());
            // printf("u2: %.6f .. %.6f, v2: %.6f .. %.6f\n", u2.min(), u2.max(), v2.min(), v2.max());
            // printf("********\n");
        }
    }
}

// biDir optical flow
void OpticalFlow::biIRLS(DImage &du, DImage &dv,
                         const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                         const DImage &u, const DImage &v,
                         const DImage &ur, const DImage &vr,
                         double as, double ap, int nIRLSIter, int nSORIter)
{
    int width = Ix.nWidth(), height = Ix.nHeight();
    DImage psi_1st, phi_1st, Ix2, Iy2, Ixt, Iyt, Ixy, lapU, lapV;
    DImage ix2, iy2, ixt, iyt, ixy, A11, A22, A12, b1, b2;
    DImage wur, wvr, urx, ury, vrx, vry, uu, vv;
    DImage wrx2(width, height), wry2(width, height);
    DImage wrxy(width, height), wrxt(width, height), wryt(width, height);
    double utmp, vtmp, gamma;
    
    du.create(width, height);//match size and set to 0
    dv.create(width, height);
    
    // symetric term
    warpImage(wur, ur, ur, u, v);
    warpImage(wvr, vr, vr, u, v);
    grad1st(urx, ury, wur);
    grad1st(vrx, vry, wvr);
        
    for (int irls = 0; irls < nIRLSIter; ++irls)
    {
        add(uu, u, du);// uu = u + du
        add(vv, v, dv);// vv = v + dv
        
        phi_d(phi_1st, uu, vv);
        psi_d(psi_1st, Ix, Iy, It, du, dv);
        
        for (int i = 0; i < u.nElements(); ++i)
        {
            utmp = u[i]+du[i]+wur[i] + urx[i]*du[i] + ury[i]*dv[i];
            vtmp = v[i]+dv[i]+wvr[i] + vrx[i]*du[i] + vry[i]*dv[i];
            gamma = ap * 0.5 / sqrt(utmp*utmp + vtmp*vtmp + 0.1);

            wrx2[i] = gamma * ((urx[i]+1)*(urx[i]+1) + (vrx[i]*vrx[i]));
            wry2[i] = gamma * ((vry[i]+1)*(vry[i]+1) + (ury[i]*ury[i]));
            wrxy[i] = gamma * ((urx[i]+1)*ury[i] + vrx[i]*(vry[i]+1));
            wrxt[i] = gamma * ((urx[i]+1)*(u[i]+wur[i]) + vrx[i]*(v[i]+wvr[i]));
            wryt[i] = gamma * ((vry[i]+1)*(v[i]+wvr[i]) + ury[i]*(u[i]+wur[i]));
        }
        
        // components for linear system
        multiply(Ix2, Ix, Ix, psi_1st);
        collapse(ix2, Ix2);
        
        multiply(Iy2, Iy, Iy, psi_1st);
        collapse(iy2, Iy2);
        
        multiply(Ixy, Ix, Iy, psi_1st);
        collapse(ixy, Ixy);
        
        multiply(Ixt, It, Ix, psi_1st);
        collapse(ixt, Ixt);
        
        multiply(Iyt, It, Iy, psi_1st);
        collapse(iyt, Iyt);

        // apply mask
        multiply(phi_1st, mask);
        
        multiply(A11, ix2, mask);
        multiply(wrx2, mask);
        add(A11, wrx2);
        add(A11, as*0.05); // add epsilon to avoid dividing zero

        multiply(A22, iy2, mask);
        multiply(wry2, mask);        
        add(A22, wry2);
        add(A22, as*0.05); // add epsilon to avoid dividing zero

        multiply(A12, ixy, mask);
        multiply(wrxy, mask);
        add(A12, wrxy);

        multiply(b1, ixt, mask);
        multiply(wrxt, mask);
        add(b1, wrxt);

        multiply(b2, iyt, mask);
        multiply(wryt, mask);
        add(b2, wryt);
        
        weighted_lap(lapU, u, phi_1st);
        weighted_lap(lapV, v, phi_1st);
        for (int i = 0; i < iyt.nElements(); ++i)
        {
            b1[i] = -b1[i] + as * lapU[i];
            b2[i] = -b2[i] + as * lapV[i];
        }
        
        // SOR iteration
        du.set(0);
        dv.set(0);
        const double omega = 1.8;
        double l, l_du, l_dv;
        int offset, tmp;

        for (int siter = 0; siter < nSORIter; siter++)
        {
            for (int h = 0; h < height; ++h)
            {
                for (int w = 0; w < width; ++w)
                {
                    offset = h * width + w;
                    l_du = 0, l_dv = 0, l = 0;
                        
                    if (h > 0)
                    {
                        tmp = offset - width;
                        l_du += phi_1st[tmp] * du[tmp];
                        l_dv += phi_1st[tmp] * dv[tmp];
                        l -= phi_1st[tmp];
                    }
                    if (h < height-1)
                    {
                        tmp = offset + width;
                        l_du += phi_1st[offset] * du[tmp];
                        l_dv += phi_1st[offset] * dv[tmp];
                        l -= phi_1st[offset];
                    }
                    if (w > 0)
                    {
                        tmp = offset - 1;
                        l_du += phi_1st[tmp] * du[tmp];
                        l_dv += phi_1st[tmp] * dv[tmp];
                        l -= phi_1st[tmp];
                    }
                    if (w < width-1)
                    {
                        tmp = offset + 1;
                        l_du += phi_1st[offset] * du[tmp];
                        l_dv += phi_1st[offset] * dv[tmp];
                        l -= phi_1st[offset];
                    }

                    l *= as;
                    l_du *= as;
                    l_dv *= as;
                        
                    // du
                    l_du = b1[offset] + l_du - A12[offset]*dv[offset];
                    du[offset] = (1-omega)*du[offset]+omega/(A11[offset]-l)*l_du;
                        
                    // dv
                    l_dv = b2[offset] + l_dv - A12[offset]*du[offset];
                    dv[offset] = (1-omega)*dv[offset]+omega/(A22[offset]-l)*l_dv;
                }
            }
        }
//        printf("SOR: du: %.6f .. %.6f, dv: %.6f .. %.6f\n", du.min(), du.max(), dv.min(), dv.max());
    }
}

// adaptive optical flow
void OpticalFlow::adIRLS(DImage &du, DImage &dv,
                         const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                         const DImage &u, const DImage &v,
                         const DImage &ur, const DImage &vr,
                         double as, double ap, int nIRLSIter, int nSORIter)
{
    int width = Ix.nWidth(), height = Ix.nHeight(), channels = Ix.nChannels();
    DImage Ix2, Iy2, Ixt, Iyt, Ixy, lapU, lapV, ix2, iy2, ixt, iyt, ixy;
    DImage wur, wvr, urx, ury, vrx, vry, uu, vv, ux, uy, vx, vy;
    DImage wrx2(width, height), wry2(width, height);
    DImage wrxy(width, height), wrxt(width, height), wryt(width, height);
    DImage psid(width, height, channels), phid(width, height);
    DImage thetad(width, height), D(width, height);
    DImage A11, A22, A12, b1, b2;
    
    du.create(width, height);//match size and set to 0
    dv.create(width, height);
    
    // symetric term
    warpImage(wur, ur, ur, u, v);
    warpImage(wvr, vr, vr, u, v);
    grad1st(urx, ury, wur);
    grad1st(vrx, vry, wvr);
        
    for (int irls = 0; irls < nIRLSIter; ++irls)
    {
        add(uu, u, du);// uu = u + du
        add(vv, v, dv);// vv = v + dv

        phi_d(phid, uu, vv);
        psi_d(psid, Ix, Iy, It, du, dv);

        // compute coefficient matrix D, ap * theta', wrx2, wry2, wrxy, wrxt, wryt
        const double kapa = 1;
        const double one = 1;
        for (int i = 0; i < u.nElements(); ++i)
        {
            double utmp = u[i]+du[i]+wur[i] + urx[i]*du[i] + ury[i]*dv[i];
            double vtmp = v[i]+dv[i]+wvr[i] + vrx[i]*du[i] + vry[i]*dv[i];
            double sym = utmp*utmp + vtmp*vtmp;

            D[i] = one / sqrt(one + kapa * sym);
            thetad[i] = ap * 0.5 / sqrt(sym + 0.1);
            
            wrx2[i] = (urx[i]+1)*(urx[i]+1) + (vrx[i]*vrx[i]);
            wry2[i] = (vry[i]+1)*(vry[i]+1) + (ury[i]*ury[i]);
            wrxy[i] = (urx[i]+1)*ury[i] + vrx[i]*(vry[i]+1);
            wrxt[i] = (urx[i]+1)*(u[i]+wur[i]) + vrx[i]*(v[i]+wvr[i]);
            wryt[i] = (vry[i]+1)*(v[i]+wvr[i]) + ury[i]*(u[i]+wur[i]);
        }
               
        // components for linear system
        // printf("D:  %.6f .. %.6f\n", D.min(), D.max());
        
        multiply(Ix2, Ix, Ix, psid);
        collapse(ix2, Ix2);
        multiply(ix2, D);
        
        multiply(Iy2, Iy, Iy, psid);
        collapse(iy2, Iy2);
        multiply(iy2, D);
        
        multiply(Ixy, Ix, Iy, psid);
        collapse(ixy, Ixy);
        multiply(ixy, D);
        
        multiply(Ixt, It, Ix, psid);
        collapse(ixt, Ixt);
        multiply(ixt, D);
        
        multiply(Iyt, It, Iy, psid);
        collapse(iyt, Iyt);
        multiply(iyt, D);

        weighted_lap(lapU, u, phid);
        weighted_lap(lapV, v, phid);
        
        multiply(wrx2, thetad);
        multiply(wrxy, thetad);
        multiply(wry2, thetad);
        multiply(wrxt, thetad);
        multiply(wryt, thetad);

        add(A11, wrx2, ix2);
        add(A22, wry2, iy2);
        add(A12, wrxy, ixy);
        add(b1, wrxt, ixt);
        substract(b1, lapU, as); // b1 -= as * lapU
        add(b2, wryt, iyt);
        substract(b2, lapV, as); // b2 -= as * lapV
        
        // apply mask
        multiply(phid, mask);
        multiply(A11, mask);
        multiply(A22, mask);
        multiply(A12, mask);
        multiply(b1, mask);
        multiply(b2, mask);
        
        add(A11, as*0.05); // add epsilon to avoid dividing zero        
        add(A22, as*0.05); // add epsilon to avoid dividing zero
        
        // SOR iteration
        du.set(0);
        dv.set(0);
        const double omega = 1.8;
        double l, l_du, l_dv;
        int offset, tmp;

        for (int siter = 0; siter < nSORIter; siter++)
        {
            for (int h = 0; h < height; ++h)
            {
                for (int w = 0; w < width; ++w)
                {
                    offset = h * width + w;
                    l_du = 0, l_dv = 0, l = 0;
                        
                    if (h > 0)
                    {
                        tmp = offset - width;
                        l_du += phid[tmp] * du[tmp];
                        l_dv += phid[tmp] * dv[tmp];
                        l -= phid[tmp];
                    }
                    if (h < height-1)
                    {
                        tmp = offset + width;
                        l_du += phid[offset] * du[tmp];
                        l_dv += phid[offset] * dv[tmp];
                        l -= phid[offset];
                    }
                    if (w > 0)
                    {
                        tmp = offset - 1;
                        l_du += phid[tmp] * du[tmp];
                        l_dv += phid[tmp] * dv[tmp];
                        l -= phid[tmp];
                    }
                    if (w < width-1)
                    {
                        tmp = offset + 1;
                        l_du += phid[offset] * du[tmp];
                        l_dv += phid[offset] * dv[tmp];
                        l -= phid[offset];
                    }

                    l *= as;
                    l_du *= as;
                    l_dv *= as;
                        
                    // du
                    l_du = -b1[offset] + l_du - A12[offset]*dv[offset];
                    du[offset] = (1-omega)*du[offset]+omega/(A11[offset]-l)*l_du;
                        
                    // dv
                    l_dv = -b2[offset] + l_dv - A12[offset]*du[offset];
                    dv[offset] = (1-omega)*dv[offset]+omega/(A22[offset]-l)*l_dv;
                }
            }
        }
        // printf("SOR: du: %.6f .. %.6f, dv: %.6f .. %.6f\n", du.min(), du.max(), dv.min(), dv.max());
    }
}

// coarse to fine strategy for optical flow
void OpticalFlow::adC2FFlow(DImage &u, DImage &v, 
                            const DImage &im1, const DImage &im2,
                            const DImage &mask1, const DImage &mask2,
                            double as, double ratio, int minWidth,
                            int nOutIter, int nIRLSIter, int nSORIter)
{
    GaussianPyramid pyr1, pyr2, mpyr1, mpyr2;

    pyr1.build(im1, ratio, minWidth);
    pyr2.build(im2, ratio, minWidth);

    // do not smooth while build pyramid
    mpyr1.build(mask1, ratio, minWidth);
    mpyr2.build(mask2, ratio, minWidth);

    // printf("pyramid constructed\n");
    const int wsize = 5;
    const double truncate = 0.1;
    DImage fIm1, fIm2, warpI2, Ix, Iy, It, mask, du, dv, tmp, D, covUV;
    int width, height;
    double xfactor, yfactor, maxVal;
    char buf[256];
    const char *outImg = "/home/iaml/Projects/exp/%s%03d.jpg";
    
    // iterate from the top level to the bottom
    for (int k = pyr1.nLevels()-1; k >= 0; --k)
    {
        // printf("Pyramid level %d\n", k);

        width = pyr1[k].nWidth();
        height = pyr1[k].nHeight();
        im2feature(fIm1, pyr1[k]);
        im2feature(fIm2, pyr2[k]);

        // if at the top level
        if (k == pyr1.nLevels()-1)
        {
            u.create(width, height);
            v.create(width, height);
            fIm2.copyTo(warpI2);
        } else {
            xfactor = (double)width / pyr1[k+1].nWidth();
            yfactor = (double)height / pyr1[k+1].nHeight();
            
            imresize(tmp, u, width, height);
            multiply(u, tmp, xfactor);
            imresize(tmp, v, width, height);
            multiply(v, tmp, yfactor);
            warpImage(warpI2, fIm1, fIm2, u, v);
        }
        
        mask.create(width, height, 1, 1);
        D.create(width, height);
        for (int l = 0; l < nOutIter+k; l++)
        {
            getGrads(Ix, Iy, It, fIm1, warpI2);
            if (l >= nOutIter+k-3)
                genInImageMask(mask, mpyr1[k], mpyr2[k], u, v);

            // co-variance orientation, magnitued
            covariance(covUV, u, v, wsize);
            maxVal = -1;
            for (int i = 0; i < covUV.nElements(); ++i)
            {
                covUV[i] = fabs(covUV[i]);
                if (maxVal < covUV[i]) maxVal = covUV[i];
            }

            // normalize
            if (maxVal < ESP)
                D.set(1);
            else
            {
                for (int i = 0; i < covUV.nElements(); ++i)
                {
                    D[i] = covUV[i] / maxVal;
                    if (D[i] <= truncate) D[i] = 0;
                    D[i] = 1 - D[i];
                }
            }
            
            adIRLS2(du, dv, D, Ix, Iy, It, mask, u, v, as, ap, nIRLSIter, nSORIter+k*3);
            
            add(u, du);
            add(v, dv);
            warpImage(warpI2, fIm1, fIm2, u, v);
            
            // for test purpose
            if (k == 0)
            {
                // show confidence map
                printf("co-variance: %.6f .. %.6f\n", covUV.min(), covUV.max());
                printf("confidence: %.6f .. %.6f\n", D.min(), D.max());

                sprintf(buf, outImg, "confidence", l);
                imwrite(buf, D);

                UCImage ucimg;
                flow2color(ucimg, u, v);
                sprintf(buf, outImg, "flow", l);
                imwrite(buf, ucimg);

                coverLabels(tmp, im1, D);
                sprintf(buf, outImg, "merge", l);
                imwrite(buf, tmp);
            }
        }
    }
}

// adaptive optical flow
void OpticalFlow::adIRLS2(DImage &du, DImage &dv, const DImage &D,
                         const DImage &Ix, const DImage &Iy, const DImage &It, const DImage &mask,
                         const DImage &u, const DImage &v,
                         double as, int nIRLSIter, int nSORIter)
{
    int width = Ix.nWidth(), height = Ix.nHeight(), channels = Ix.nChannels();
    DImage Ix2, Iy2, Ixt, Iyt, Ixy, lapU, lapV, ix2, iy2, ixt, iyt, ixy, uu, vv;
    DImage psid(width, height, channels), phid(width, height);
    DImage A11, A22, A12, b1, b2;
    
    du.create(width, height);//match size and set to 0
    dv.create(width, height);
    
    for (int irls = 0; irls < nIRLSIter; ++irls)
    {
        add(uu, u, du);// uu = u + du
        add(vv, v, dv);// vv = v + dv

        phi_d(phid, uu, vv);
        multiply(phid, D);
        psi_d(psid, Ix, Iy, It, du, dv);

        multiply(Ix2, Ix, Ix, psid);
        collapse(ix2, Ix2);
        
        multiply(Iy2, Iy, Iy, psid);
        collapse(iy2, Iy2);
        
        multiply(Ixy, Ix, Iy, psid);
        collapse(ixy, Ixy);
        
        multiply(Ixt, It, Ix, psid);
        collapse(ixt, Ixt);
        
        multiply(Iyt, It, Iy, psid);
        collapse(iyt, Iyt);

        weighted_lap(lapU, u, phid);
        weighted_lap(lapV, v, phid);

        multiply(A11, ix2, D);
        multiply(A22, iy2, D);
        multiply(A12, ixy, D);
        multiply(b1, ixt, D);
        substract(b1, lapU, as); // b1 -= as * lapU
        multiply(b2, iyt, D);
        substract(b2, lapV, as); // b2 -= as * lapV
        
        // apply mask
        multiply(phid, mask);
        multiply(A11, mask);
        multiply(A22, mask);
        multiply(A12, mask);
        multiply(b1, mask);
        multiply(b2, mask);
        
        add(A11, as*0.05); // add epsilon to avoid dividing zero        
        add(A22, as*0.05); // add epsilon to avoid dividing zero
        
        // SOR iteration
        du.set(0);
        dv.set(0);
        const double omega = 1.8;
        double l, l_du, l_dv;
        int offset, tmp;

        for (int siter = 0; siter < nSORIter; siter++)
        {
            for (int h = 0; h < height; ++h)
            {
                for (int w = 0; w < width; ++w)
                {
                    offset = h * width + w;
                    l_du = 0, l_dv = 0, l = 0;
                        
                    if (h > 0)
                    {
                        tmp = offset - width;
                        l_du += phid[tmp] * du[tmp];
                        l_dv += phid[tmp] * dv[tmp];
                        l -= phid[tmp];
                    }
                    if (h < height-1)
                    {
                        tmp = offset + width;
                        l_du += phid[offset] * du[tmp];
                        l_dv += phid[offset] * dv[tmp];
                        l -= phid[offset];
                    }
                    if (w > 0)
                    {
                        tmp = offset - 1;
                        l_du += phid[tmp] * du[tmp];
                        l_dv += phid[tmp] * dv[tmp];
                        l -= phid[tmp];
                    }
                    if (w < width-1)
                    {
                        tmp = offset + 1;
                        l_du += phid[offset] * du[tmp];
                        l_dv += phid[offset] * dv[tmp];
                        l -= phid[offset];
                    }

                    l *= as;
                    l_du *= as;
                    l_dv *= as;
                        
                    // du
                    l_du = -b1[offset] + l_du - A12[offset]*dv[offset];
                    du[offset] = (1-omega)*du[offset]+omega/(A11[offset]-l)*l_du;
                        
                    // dv
                    l_dv = -b2[offset] + l_dv - A12[offset]*du[offset];
                    dv[offset] = (1-omega)*dv[offset]+omega/(A22[offset]-l)*l_dv;
                }
            }
        }
    }
}

// coarse to fine strategy for optical flow
void OpticalFlow::stC2FFlow(std::vector<DImage> &u, std::vector<DImage> &v, 
                            const std::vector<DImage> &im,
                            const std::vector<DImage> &mask, int i0,
                            double as, double at, double ratio, int minWidth,
                            int nOutIter, int nIRLSIter, int nSORIter)
{
    assert(im.size() == mask.size() && im.size() == 3 && u.size() == v.size() &&
           v.size() == 2);

    int i1 = (i0 + 1) % 3;
    int i2 = (i1 + 1) % 3;
    DImage fIm1, fIm2, warpI2, du, dv;

    // im0 -> im1
    if (u[i0].isEmpty() || v[i0].isEmpty())
    {
        printf("previous flow is not available\n");
        adC2FFlow(u[i0], v[i0], im[i0], im[i1], mask[i0], mask[i1],
                  as, ratio, minWidth, nOutIter, nIRLSIter, nSORIter);
    }

    // im1 -> im2
    adC2FFlow(u[i1], v[i1], im[i1], im[i2], mask[i1], mask[i2],
              as, ratio, minWidth, nOutIter, nIRLSIter, nSORIter);

    // temporal smooth im0 -> im1
    temporalSmooth(u[i0], v[i0], u[i1], v[i1], im[i0], im[i1], mask[i0], mask[i1],
                   as, at, nOutIter, nIRLSIter, nSORIter);
}

void OpticalFlow::temporalSmooth(DImage &u0, DImage &v0,
                                 const DImage &u1, const DImage &v1,
                                 const DImage &im0, const DImage &im1,
                                 const DImage &mask0, const DImage &mask1,
                                 double as, double at,
                                 int nOutIter, int nIRLSIter, int nSORIter)
{
    assert(im0.match3D(im1) && mask0.match3D(mask1) && im0.match3D(mask0) &&
           u1.match3D(u) && v1.match3D(v0) && u0.match3D(v0) && u0.match2D(im0));

    int width = im0.nWidth(), height = im0.nHeight();
    
    DImage fIm0, fIm1, mask, warp;
    DImage Ix, Iy, It, Ix2, Iy2, Ixy, Ixt, Iyt, ix2, iy2, ixt, iyt, ixy;
    DImage psid, phid, thetad, lapU, lapV;
    DImage uu, vv, warpu, warpv, du(width, height), dv(width, height);
    DImage A11, A12, A22, b1, b2;
    
    im2feature(fIm0, im0);
    im2feature(fIm1, im1);
    warpImage(warp, fIm0, fIm1, u0, v0);

    for (int oiter = 0; oiter < nOutIter; ++oiter)
    {
        getGrads(Ix, Iy, It, fIm0, warp);
        du.set(0), dv.set(0);

        for (int irlsiter = 0; irlsiter < nIRLSIter; ++irlsiter)
        {
            add(uu, u, du);
            add(vv, v, dv);
            warpImage(warpu, u0, u1, u0, v0);
            warpImage(warpv, v0, v1, u0, v0);

            phi_d(phid, uu, vv);
            psi_d(psid, Ix, Iy, It, du, dv);
            
            // components for linear system
            multiply(Ix2, Ix, Ix, psid);
            collapse(ix2, Ix2);

            multiply(Iy2, Iy, Iy, psid);
            collapse(iy2, Iy2);
            
            multiply(Ixy, Ix, Iy, psid);
            collapse(ixy, Ixy);
            
            multiply(Ixt, It, Ix, psid);
            collapse(ixt, Ixt);
            
            multiply(Iyt, It, Iy, psid);
            collapse(iyt, Iyt);
            
            weighted_lap(lapU, u, phid);
            weighted_lap(lapV, v, phid);

            for (int i = 0; i < iyt.nElements(); ++i)
            {
                ixt[i] = -ixt[i] + as * lapU[i];
                iyt[i] = -iyt[i] + as * lapV[i];
            }
            
            // SOR iteration
            du.set(0);
            dv.set(0);
            const double omega = 1.8;
            double l, l_du, l_dv;
            int offset, tmp;

            for (int siter = 0; siter < nSORIter; siter++)
            {
                for (int h = 0; h < height; ++h)
                {
                    for (int w = 0; w < width; ++w)
                    {
                        offset = h * width + w;
                        l_du = 0, l_dv = 0, l = 0;
                        
                        if (h > 0)
                        {
                            tmp = offset - width;
                            l_du += phi_1st[tmp] * du[tmp];
                            l_dv += phi_1st[tmp] * dv[tmp];
                            l -= phi_1st[tmp];
                        }
                        if (h < height-1)
                        {
                            tmp = offset + width;
                            l_du += phi_1st[offset] * du[tmp];
                            l_dv += phi_1st[offset] * dv[tmp];
                            l -= phi_1st[offset];
                        }
                        if (w > 0)
                        {
                            tmp = offset - 1;
                            l_du += phi_1st[tmp] * du[tmp];
                            l_dv += phi_1st[tmp] * dv[tmp];
                            l -= phi_1st[tmp];
                        }
                        if (w < width-1)
                        {
                            tmp = offset + 1;
                            l_du += phi_1st[offset] * du[tmp];
                            l_dv += phi_1st[offset] * dv[tmp];
                            l -= phi_1st[offset];
                        }

                        l *= as;
                        l_du *= as;
                        l_dv *= as;
                        
                        // du
                        l_du = ixt[offset] + l_du - ixy[offset] * dv[offset];
                        du[offset] = (1-omega)*du[offset]+omega/(ix2[offset]-l)*l_du;
                        
                        // dv
                        l_dv = iyt[offset] + l_dv - ixy[offset] * du[offset];
                        dv[offset] = (1-omega)*dv[offset]+omega/(iy2[offset]-l)*l_dv;
                    }
                }

            }
//            printf("du: %.6f .. %.6f, dv: %.6f .. %.6f\n", du.min(), du.max(), dv.min(), dv.max());
        }

        add(u, du);// u += du
        add(v, dv);// v += dv
        warpImage(warp, im1, im2, u, v);
    }
            
        }
    }
}
