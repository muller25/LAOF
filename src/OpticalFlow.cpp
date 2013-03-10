#include "OpticalFlow.h"

void OpticalFlow::warpImage(DImage &warp, const DImage &im1, const DImage &im2,
                            const DImage &u, const DImage &v)
{
    assert(im1.match3D(im2) && u.match3D(v) && im1.match2D(u));
    
    int width = im1.nWidth(), height = im1.nHeight(), channels = im1.nChannels();
    double nx, ny;
    int offset;

    warp.create(width, height, channels);    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            nx = w + u[offset];
            ny = h + v[offset];
            offset *= channels;
            
            if (nx < 0 || nx > width-1 || ny < 0 || ny > height-1)
            {
                for (int k = 0; k < channels; k++)
                    warp[offset+k] = im1[offset+k];

                continue;
            }

            biInterpolate(warp.ptr()+offset, im2, nx, ny);
        }
    }
}

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
                if (lapPara[k] < 1E-20) continue;

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
    double filter[5] = {0.02, 0.11, 0.74, 0.11, 0.02};

    DImage sIm1, sIm2, im;
    filtering(sIm1, im1, filter, 2, filter, 2);
    filtering(sIm2, im2, filter, 2, filter, 2);

    addWeighted(im, sIm1, 0.4, sIm2, 0.6);
    grad1st(Ix, Iy, im);
    substract(It, sIm2, sIm1);
}

void OpticalFlow::estLapNoise(const DImage &im1, const DImage &im2)
{
    assert(im1.match3D(im2));

    int height = im1.nHeight(), width = im1.nWidth(), channels = im1.nChannels();
    int offset;
    std::vector<double> total(channels, 0);
    double tmp;
    
    lapPara.assign(channels, 0);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            for (int k = 0; k < channels; ++k)
            {
                tmp = fabs(im1[offset+k] - im2[offset+k]);
                if (tmp > 0 && tmp < 1000000)
                {
                    lapPara[k] += tmp;
                    total[k]++;
                }
            }
        }
    }
    
    for (int k = 0; k < channels; ++k)
    {
        if (total[k] == 0)
        {
            printf("All the pixels are invalid in estimation Laplacian noise!!!\n");
            printf("Something severely wrong happened!!!\n");
			lapPara[k] = 0.001;
		}
		else
			lapPara[k] /= total[k];
    }
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

    // if (channels == 3)
    // {
    //     gradX(gx, im);
    //     gradY(gy, im);
    //     pgx = gx.ptr(), pgy = gy.ptr();
        
    //     // mix channels
    //     for (int h = 0; h < height; ++h)
    //     {
    //         for (int w = 0; w < width; ++w)
    //         {
    //             offset = h * width + w;
    //             foffset = offset * nchannels;
    //             offset *= channels;
    //             for (int k = 0; k < channels; ++k)
    //             {
    //                 pf[foffset + k] = pi[offset + k];
    //                 pf[foffset + k + channels] = gamma * pgx[offset + k];
    //                 pf[foffset + k + 2*channels] = gamma * pgx[offset + k];
    //             }
    //         }
    //     }
        
    //     return;
    // }

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
                pf[foffset + 3] = im[offset*3+1] - im[offset*3];
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
            x = nu, y = nv;
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
    
    // init lap noise
    lapPara.assign(im1.nChannels()+2, 0.02);

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
            imresize(tmp, u, width, height);
            multiply(u, tmp, 1./ratio);

            imresize(tmp, v, width, height);
            multiply(v, tmp, 1./ratio);

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

        // esitmate noise level
        estLapNoise(im1, warp);
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
    mpyr1.build(mask1, ratio, minWidth);
    mpyr2.build(mask2, ratio, minWidth);
    printf("pyramid constructed\n");
    
    DImage fIm1, fIm2, warpI1, warpI2, Ix, Iy, It, mask;
    DImage du1, dv1, du2, dv2, tmp;
    int width, height;
    
    // init lap noise
    lapPara.assign(im1.nChannels()+2, 0.02);

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
            u1.create(width, height);
            v1.create(width, height);
            fIm2.copyTo(warpI2);
            
            u2.create(width, height);
            v2.create(width, height);
            fIm1.copyTo(warpI1);
        } else {
            imresize(tmp, u1, width, height);
            multiply(u1, tmp, 1./ratio);
            imresize(tmp, v1, width, height);
            multiply(v1, tmp, 1./ratio);
            warpImage(warpI2, fIm1, fIm2, u1, v1);
            
            imresize(tmp, u2, width, height);
            multiply(u2, tmp, 1./ratio);
            imresize(tmp, v2, width, height);
            multiply(v2, tmp, 1./ratio);
            warpImage(warpI1, fIm2, fIm1, u2, v2);
        }
        
        for (int l = 0; l < nBiIter; l++)
        {
            // forward flow
            getGrads(Ix, Iy, It, fIm1, warpI2);
            genInImageMask(mask, mpyr1[k], mpyr2[k], u1, v1);
//            biIRLS(du1, dv1, Ix, Iy, It, mask, u1, v1, u2, v2, as, ap, nIRLSIter, nSORIter);
            adIRLS(du1, dv1, Ix, Iy, It, mask, u1, v1, u2, v2, as, ap, nIRLSIter, nSORIter);

            // backward flow
            getGrads(Ix, Iy, It, fIm2, warpI1);
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
//        multiply(ix2, mask);
//        gaussianBlur(A11, ix2, 1.2, 3);
        multiply(wrx2, mask);
        add(A11, wrx2);
        add(A11, as*0.05); // add epsilon to avoid dividing zero

        multiply(A22, iy2, mask);
//        multiply(iy2, mask);
//        gSmooth(A22, iy2, 1.2, 3);
        multiply(wry2, mask);        
        add(A22, wry2);
        add(A22, as*0.05); // add epsilon to avoid dividing zero

        multiply(A12, ixy, mask);
//        multiply(ixy, mask);
//        gSmooth(A12, ixy, 1.2, 3);
        multiply(wrxy, mask);
        add(A12, wrxy);

        multiply(b1, ixt, mask);
//        multiply(ixt, mask);
//        gSmooth(b1, ixt, 1.2, 3);
        multiply(wrxt, mask);
        add(b1, wrxt);

        multiply(b2, iyt, mask);
//        multiply(iyt, mask);
//        gSmooth(b2, iyt, 1.2, 3);
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

// spatio-temporal optical flow
// given image im[t-1], im[t], im[t+1], im[t+2], flow f
// f[t-1]: im[t-1] -> im[t], f[t]: im[t] -> im[t+1], f[t+1]: im[t+1] -> im[t+2]
// reverse flow rf
// f[t-1]: im[t] -> im[t-1], f[t]: im[t+1] -> im[t], f[t+1]: im[t+2] -> im[t+1]
void OpticalFlow::stFlow(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
                         const DImage &im1, const DImage &im2,
                         const DImage &mask1, const DImage &mask2,
                         const DImage &pu, const DImage &pv,
                         const DImage &nu, const DImage &nv,
                         const DImage &pur, const DImage &pvr,
                         const DImage &nur, const DImage &nvr,
                         double as, double ap, int nBiIter, int nIRLSIter, int nSORIter)
{
    DImage Ix, Iy, It, du1, dv1, du2, dv2, pphid, pphidr, warpI1, warpI2, mask;
    DImage fIm1, fIm2, wpu, wpv, wpur, wpvr, wnu, wnv, tmpu, tmpv;

    im2feature(fIm1, im1);
    im2feature(fIm2, im2);

    fIm1.copyTo(warpI1);
    fIm2.copyTo(warpI2);

    // u(p-w) = u[t-1]
    multiply(tmpu, pu, -1);
    multiply(tmpv, pv, -1);
    warpImage(wpu, pu, pu, tmpu, tmpv);
    warpImage(wpv, pv, pv, tmpu, tmpv);
    phi_d(pphid, wpu, wpv);

    // ur(p-w) = ur[t-1]
    multiply(tmpu, pur, -1);
    multiply(tmpv, pvr, -1);
    warpImage(wpur, pur, pur, tmpu, tmpv);
    warpImage(wpvr, pvr, pvr, tmpu, tmpv);
    phi_d(pphidr, wpur, wpvr);
    for (int i = 0; i < nBiIter; ++i)
    {
        // forward flow: im1 -> im2
        warpImage(wnu, nu, nu, u1, v1);
        warpImage(wnv, nv, nv, u1, v1);
    
        getGrads(Ix, Iy, It, fIm1, warpI2);
        genInImageMask(mask, mask1, mask2, u1, v1);
        adIRLS3(du1, dv1, Ix, Iy, It, mask, pphid, wpu, wpv, u1, v1, wnu, wnv, u2, v2,
                as, ap, nIRLSIter, nSORIter);

        // backward flow: cmask -> pmask
        warpImage(wnu, nur, nur, u2, v2);
        warpImage(wnv, nvr, nvr, u2, v2);
        
        getGrads(Ix, Iy, It, fIm2, warpI1);
        genInImageMask(mask, mask2, mask1, u2, v2);
        adIRLS3(du2, dv2, Ix, Iy, It, mask, pphidr, wpur, wpvr, u2, v2, wnu, wnv, u1, v1,
                as, ap, nIRLSIter, nSORIter);

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

// adaptive, spatio-temporal optical flow
void OpticalFlow::adIRLS3(DImage &du, DImage &dv,
                          const DImage &Ix, const DImage &Iy, const DImage &It,
                          const DImage &mask, const DImage &pphid,
                          const DImage &pu, const DImage &pv,
                          const DImage &u, const DImage &v,
                          const DImage &nu, const DImage &nv,
                          const DImage &ur, const DImage &vr,
                          double as, double ap, int nIRLSIter, int nSORIter)
{
    int width = Ix.nWidth(), height = Ix.nHeight();
    DImage Ix2, Iy2, Ixt, Iyt, Ixy, lapU, lapV, ix2, iy2, ixt, iyt, ixy, psid, phid;
    DImage wur, wvr, urx, ury, vrx, vry, uu, vv, ux, uy, vx, vy;
    DImage wrx2(width, height), wry2(width, height);
    DImage wrxy(width, height), wrxt(width, height), wryt(width, height);
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

        // spatio-temporal
        weighted_lap3(lapU, pu, u, nu, pphid, phid);

        // spatio-temporal
        weighted_lap3(lapV, pu, v, nu, pphid, phid); 
        
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

                    // spatio-temporal divergence
                    l -= (phid[offset] + pphid[offset] * mask[offset]);
                    
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
