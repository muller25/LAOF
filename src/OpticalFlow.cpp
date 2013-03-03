#include "OpticalFlow.h"

void OpticalFlow::warpImage(DImage &warp, const DImage &im1, const DImage &im2,
                            const DImage &u, const DImage &v)
{
    assert(im1.match3D(im2) && u.match3D(v) && im1.match2D(u));
    
    int width = im1.nWidth(), height = im1.nHeight(), channels = im1.nChannels();
    double *p1 = im1.ptr();
    
    warp.create(width, height, channels);
    double *pw = warp.ptr();

    double *pu = u.ptr();
    double *pv = v.ptr();
    double nx, ny;
    int offset;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            nx = w + pu[offset];
            ny = h + pv[offset];
            offset *= channels;
            
            if (nx < 0 || nx > width-1 || ny < 0 || ny > height-1)
            {
                for (int k = 0; k < channels; k++)
                    pw[offset+k] = p1[offset+k];

                continue;
            }

            biInterpolate(pw+offset, im2, nx, ny);
        }
    }
}

// 计算im1到im2的光流，u, v为初始化值，并将结果保存在u,v。a_s是smoothness term的权重
// u, v, warp需要预先初始化
void OpticalFlow::SORSolver(DImage &u, DImage &v, DImage &warp,
                            const DImage &im1, const DImage &im2,
                            double a_s, int nOutIter, int nInIter, int nSORIter)
{
    assert(im1.match3D(im2) && u.match3D(v) && warp.match3D(im1));

    int width = im1.nWidth(), height = im1.nHeight();
    DImage du(width, height), dv(width, height);
    double *pdu = du.ptr(), *pdv = dv.ptr();
    DImage It, Ix, Iy, Ixy, Ix2, Iy2, Ixt, Iyt, uu, vv;
    DImage psi_1st, phi_1st, lapU, lapV;
    DImage ixy, ix2, iy2, iyt, ixt;
    double *pixy, *pix2, *piy2, *pixt, *piyt, *pphi;
    
    // outer fixed point iteration for u, v
    for (int oiter = 0; oiter < nOutIter; oiter++)
    {
//        printf("outer fixed point iterations %d\n", oiter);

        // compute Ix, Iy and It
        getGrads(Ix, Iy, It, im1, warp);

        // inner fixed point iteration for du, dv
        du.setTo(0);
        dv.setTo(0);
        for (int iiter = 0; iiter < nInIter; iiter++)
        {
            add(uu, u, du);// uu = u + du
            add(vv, v, dv);// vv = v + dv
            phi_d(phi_1st, uu, vv);
            pphi = phi_1st.ptr();
            psi_d(psi_1st, Ix, Iy, It, du, dv);
            
            // components for linear system
            multiply(Ix2, Ix, Ix, psi_1st);
            collapse(ix2, Ix2);
            pix2 = ix2.ptr();

            multiply(Iy2, Iy, Iy, psi_1st);
            collapse(iy2, Iy2);
            piy2 = iy2.ptr();
            
            multiply(Ixy, Ix, Iy, psi_1st);
            collapse(ixy, Ixy);
            pixy = ixy.ptr();
            
            multiply(Ixt, It, Ix, psi_1st);
            collapse(ixt, Ixt);
            pixt = ixt.ptr();
            
            multiply(Iyt, It, Iy, psi_1st);
            collapse(iyt, Iyt);
            piyt = iyt.ptr();
            
            weighted_lap(lapU, u, phi_1st);
            weighted_lap(lapV, v, phi_1st);
            for (int i = 0; i < iyt.nElements(); ++i)
            {
                ixt[i] = -ixt[i] - a_s * lapU[i];
                iyt[i] = -iyt[i] - a_s * lapV[i];
            }

            // SOR iteration
            du.setTo(0);
            dv.setTo(0);
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
                            l_du += pphi[tmp] * pdu[tmp];
                            l_dv += pphi[tmp] * pdv[tmp];
                            l += pphi[tmp];
                        }
                        if (h < height-1)
                        {
                            tmp = offset + width;
                            l_du += pphi[offset] * pdu[tmp];
                            l_dv += pphi[offset] * pdv[tmp];
                            l += pphi[offset];
                        }
                        if (w > 0)
                        {
                            tmp = offset - 1;
                            l_du += pphi[tmp] * pdu[tmp];
                            l_dv += pphi[tmp] * pdv[tmp];
                            l += pphi[tmp];
                        }
                        if (w < width-1)
                        {
                            tmp = offset + 1;
                            l_du += pphi[offset] * pdu[tmp];
                            l_dv += pphi[offset] * pdv[tmp];
                            l += pphi[offset];
                        }

                        l *= a_s;
                        l_du *= -a_s;
                        l_dv *= -a_s;
                        
                        // du
                        l_du = pixt[offset] - l_du - pixy[offset] * pdv[offset];
                        pdu[offset] = (1-omega)*pdu[offset]+omega/(pix2[offset]+l+a_s*0.05)*l_du;
                        
                        // dv
                        l_dv = piyt[offset] - l_dv - pixy[offset] * pdu[offset];
                        pdv[offset] = (1-omega)*pdv[offset]+omega/(piy2[offset]+l+a_s*0.05)*l_dv;
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

// calculate res = psi_1st((It + Ix * du + Iy * dv)^2)
void OpticalFlow::psi_d(DImage &res, const DImage &Ix, const DImage &Iy,
                        const DImage &It, const DImage &du, const DImage &dv)
{
    assert(Ix.match3D(Iy) && Iy.match3D(It) && du.match3D(dv) &&
           It.match2D(dv) && du.nChannels() == 1);

    const double epsilon = 1e-6;
    
    int width = Ix.nWidth(), height = Ix.nHeight(), channels = Ix.nChannels();
    double *px = Ix.ptr(), *py = Iy.ptr(), *pt = It.ptr();

    res.create(width, height, channels);
    double *p = res.ptr(), *pu = du.ptr(), *pv = dv.ptr(), tmp;
    int io, fo;
    
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            fo = h * width + w;
            io = fo * channels;
            for (int k = 0; k < channels; ++k)
            {
                if (lapPara[k] < 1E-20) continue;

                tmp = pt[io+k] + px[io+k]*pu[fo] + py[io+k]*pv[fo];
                tmp *= tmp;
                p[io+k] = 0.5 / sqrt(tmp + epsilon);
            }
        }
    }
}

// calculate res = phi_1st(u^2 + v^2)
void OpticalFlow::phi_d(DImage &res, const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const double epsilon = 1e-6;
    double tmp;
    DImage ux, uy, vx, vy;
    double *pux, *puy, *pvx, *pvy, *p;
    
    gradX(ux, u);
    gradY(uy, u);
    gradX(vx, v);
    gradY(vy, v);
    pux = ux.ptr(), puy = uy.ptr(), pvx = vx.ptr(), pvy = vy.ptr();
    
    res.create(u.nWidth(), u.nHeight());
    p = res.ptr();
    
    for (int i = 0; i < res.nElements(); i++)
    {
        tmp = pux[i]*pux[i] + puy[i]*puy[i] + pvx[i]*pvx[i] + pvy[i]*pvy[i];
        p[i] = 0.5 / sqrt(tmp + epsilon);
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
    gradX(Ix, im);
    gradY(Iy, im);
    substract(It, sIm2, sIm1);
}

void OpticalFlow::estLapNoise(const DImage &im1, const DImage &im2)
{
    assert(im1.match3D(im2));

    int height = im1.nHeight(), width = im1.nWidth(), channels = im1.nChannels();
    double *p1 = im1.ptr(), *p2 = im2.ptr(), tmp;
    int offset;
    std::vector<double> total(channels, 0);
    
    lapPara.assign(channels, 0);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            for (int k = 0; k < channels; ++k)
            {
                tmp = fabs(p1[offset+k] - p2[offset+k]);
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

void OpticalFlow::c2fFlow(DImage &u, DImage &v, const DImage &im1, const DImage &im2,
                          double a_s, double ratio, int minWidth,
                          int nOutIter, int nInIter, int nSORIter)
{
    GaussianPyramid pyr1, pyr2;

    pyr1.ConstructPyramid(im1, ratio, minWidth);
    pyr2.ConstructPyramid(im2, ratio, minWidth);

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

        SORSolver(u, v, warp, fIm1, fIm2, a_s, nOutIter+k, nInIter, nSORIter+k*3);
        printf("u: %.6f .. %.6f, v: %.6f .. %.6f\n", u.min(), u.max(), v.min(), v.max());
    }
}

void OpticalFlow::im2feature(DImage &feature, const DImage &im)
{
    int height = im.nHeight(), width = im.nWidth(), channels = im.nChannels();
    const int nchannels = channels + 2;
    DImage gx, gy;    
    double *pgx, *pgy;
    
    feature.create(width, height, nchannels);
    double *pf = feature.ptr(), *pi = im.ptr();
    int offset;
    
    if (channels == 1)
    {
        gradX(gx, im);
        gradY(gy, im);
        pgx = gx.ptr(), pgy = gy.ptr();
        
        // mix channels
        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; w++)
            {
                offset = h * width + w;
                pf[offset*nchannels] = pi[offset];
                pf[offset*nchannels+1] = pgx[offset];
                pf[offset*nchannels+2] = pgy[offset];
            }
        }
        
        return;
    }

    if (channels == 3)
    {
        DImage gray;

        desuarate(gray, im);
        double *pg = gray.ptr();
        
        gradX(gx, gray);
        gradY(gy, gray);
        pgx = gx.ptr(), pgy = gy.ptr();
        
        // mix channels
        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; ++w)
            {
                offset = h * width + w;
                pf[offset*nchannels] = pg[offset];
                pf[offset*nchannels+1] = pgx[offset];
                pf[offset*nchannels+2] = pgy[offset];
                pf[offset*nchannels+3] = pi[offset*3+1] - pi[offset*3];
                pf[offset*nchannels+4] = pi[offset*3+1] - pi[offset*3+2];
            }
        }
        
        return;
    }
    
    im.copyTo(feature);
}
