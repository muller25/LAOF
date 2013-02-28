#include "OpticalFlow.h"

/*
  [description]
  根据u，v变换im2。对于运动超出边界的像素，用im1相对的像素替换,并将结果保存到warp中。
  其中im1，im2，warp的元素应该为double类型，u,v应该为单通道浮点数。
  [params]
  im1 - 用于填充运动超出边界的像素 (in)
  im2 - 用于变换 (in)
  u - 水平运动位移 (in)
  v - 垂直运动位移 (in)
  warp - 变换后的结果 (out)
  [return]
  无
*/
void OpticalFlow::warpImage(const Mat &im1, const Mat &im2, const Mat &u, const Mat &v, Mat &warp)
{
    assert(matchAll(im1, im2) && matchAll(im2, warp) && matchAll(u, v) &&
           im1.depth() == CV_64F && u.type() == CV_64F);
    
    int rows = im1.rows, cols = im1.cols, channels = im1.channels();
    int istep = im1.step / sizeof(double);
    double *pim1 = (double *)im1.data, *pwarp = (double *)warp.data;
    
    int fstep = u.step / sizeof(double);
    double *pu = (double *)u.data, *pv = (double *)v.data;

    int r, c, k, ioffset, foffset;
    double nx, ny;

    warp.setTo(0);
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            ioffset = r * istep + c * channels;
            foffset = r * fstep + c;
            nx = c + pu[foffset];
            ny = r + pv[foffset];
                
            if (nx < 0 || nx > cols-1 || ny < 0 || ny > rows-1)
            {
                for (k = 0; k < channels; k++)
                    pwarp[ioffset+k] = pim1[ioffset+k];

                continue;
            }

            biInterpolate(im2, nx, ny, pwarp+ioffset);
        }
    }
}

// 计算im1到im2的光流，u, v为初始化值，并将结果保存在u,v。a_s是smoothness term的权重
// u, v需要预先初始化
void OpticalFlow::compute(const Mat &im1, const Mat &im2, Mat &warp,
                          Mat &u, Mat &v, const double a_s,
                          const int nOutIter, const int nInIter, const int nSORIter)
{
    int rows = im1.rows;
    int cols = im1.cols;
    int channels = im1.channels();

    Mat du(rows, cols, CV_64F), dv(rows, cols, CV_64F);
    Mat It, Ix, Iy, Ixy, Ix2, Iy2, Ixt, Iyt, uu, vv;
     
    Mat ixy(rows, cols, CV_64F), ix2(rows, cols, CV_64F);
    Mat iy2(rows, cols, CV_64F), iyt(rows, cols, CV_64F);
    Mat ixt(rows, cols, CV_64F);

    Mat psi_1st(rows, cols, CV_64FC(channels));
    Mat phi_1st(rows, cols, CV_64F);
    Mat lapU(rows, cols, CV_64F), lapV(rows, cols, CV_64F);
    
    double *pphi = (double *)phi_1st.data;
    double *pdu = (double *)du.data, *pdv = (double *)dv.data;
    double *pixy = (double *)ixy.data, *pix2 = (double *)ix2.data;
    double *piy2 = (double *)iy2.data, *pixt = (double *)ixt.data;
    double *piyt = (double *)iyt.data;
    int step = du.step / sizeof(double), offset, tmp;

    double maxu, minu;
    
    // outer fixed point iteration for u, v
    estLapNoise(im1, warp);
    for (int oiter = 0; oiter < nOutIter; oiter++)
    {
        printf("outer fixed point iterations %d\n", oiter);

        // compute Ix, Iy and It
        getGrads(im1, warp, Ix, Iy, It);
        sanityCheck(Ix, Iy, It, u, v);

        minMaxIdx(u, &minu, &maxu);
        printf("u: %.3f .. %.3f\n", minu, maxu);
        minMaxIdx(v, &minu, &maxu);
        printf("v: %.3f .. %.3f\n", minu, maxu);
        // minMaxIdx(It, &minu, &maxu);
        // printf("It: %.3f .. %.3f\n", minu, maxu);

        // inner fixed point iteration for du, dv
        du.setTo(0), dv.setTo(0);
        for (int iiter = 0; iiter < nInIter; iiter++)
        {
//            printf("inner fixed point iterations %d\n", iiter);

            add(u, du, uu);
            add(v, dv, vv);
            phi_d(uu, vv, phi_1st);
            psi_d(Ix, Iy, It, du, dv, psi_1st);

            // printf("inner fixed point iteration\n");
            // sanityCheck(Ix, Iy, It, uu, vv);

            // components for linear system
            multiply(Ix, Ix, Ix2);
            Ix2.mul(psi_1st);
            collapse<double>(Ix2, ix2);

            multiply(Iy, Iy, Iy2);
            Iy2.mul(psi_1st);
            collapse<double>(Iy2, iy2);

            multiply(Ix, Iy, Ixy);
            Ixy.mul(psi_1st);
            collapse<double>(Ixy, ixy);

            multiply(Ix, It, Ixt);
            Ixt.mul(psi_1st);
            collapse<double>(Ixt, ixt);

            multiply(Iy, It, Iyt);
            Iyt.mul(psi_1st);
            collapse<double>(Iyt, iyt);

            weighted_lap(u, phi_1st, lapU);
            addWeighted(lapU, -a_s, ixt, -1, 0, ixt);
            
            weighted_lap(v, phi_1st, lapV);
            addWeighted(lapV, -a_s, iyt, -1, 0, iyt);

            // SOR iteration
            du.setTo(0), dv.setTo(0);
            const double omega = 1.8;
            double l, l_du, l_dv;

            for (int siter = 0; siter < nSORIter; siter++)
            {
//                printf("SOR iterations %d\n", siter);
                
                for (int r = 0; r < rows; r++)
                {
                    for (int c = 0; c < cols; c++)
                    {
                        offset = r * step + c;
                        l_du = 0, l_dv = 0, l = 0;
                        
                        if (r > 0)
                        {
                            tmp = offset - step;
                            l_du += pphi[tmp] * pdu[tmp];
                            l_dv += pphi[tmp] * pdv[tmp];
                            l += pphi[tmp];
                        }
                        if (r < rows-1)
                        {
                            tmp = offset + step;
                            l_du += pphi[offset] * pdu[tmp];
                            l_dv += pphi[offset] * pdv[tmp];
                            l += pphi[offset];
                        }
                        if (c > 0)
                        {
                            tmp = offset - 1;
                            l_du += pphi[tmp] * pdu[tmp];
                            l_dv += pphi[tmp] * pdv[tmp];
                            l += pphi[tmp];
                        }
                        if (c < cols-1)
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
                        pdu[offset] = (1-omega)*pdu[offset]+omega/(pix2[offset]+l)*l_du;
                        
                        // dv
                        l_dv = piyt[offset] - l_dv - pixy[offset] * pdu[offset];
                        pdv[offset] = (1-omega)*pdv[offset]+omega/(piy2[offset]+l)*l_dv;
                    }
                }
            }
            // minMaxIdx(du, &minu, &maxu);
            // printf("du: %.3f .. %.3f\n", minu, maxu);
            // minMaxIdx(dv, &minu, &maxu);
            // printf("dv: %.3f .. %.3f\n", minu, maxu);
        }

        u += du;
        v += dv;
        warpImage(im1, im2, u, v, warp);

        // esitmate noise level
        estLapNoise(im1, warp);
    }
}

// calculate res = psi_1st((It + Ix * du + Iy * dv)^2)
void OpticalFlow::psi_d(const Mat &Ix, const Mat &Iy, const Mat &It,
                        const Mat &du, const Mat &dv, Mat &res)
{
    assert(matchAll(Ix, Iy) && matchAll(Iy, It) && matchAll(It, res) &&
           matchAll(du, dv) && Ix.depth() == CV_64F && du.depth() == CV_64F);
    
    int rows = Ix.rows, cols = Ix.cols, channels = Ix.channels();
    int istep = Ix.step / sizeof(double), io, fo;
    double *px = (double *)Ix.data, *py = (double *)Iy.data, *pt = (double *)It.data;
    double *p = (double *)res.data;
    
    int fstep = du.step / sizeof(double);
    double *pu = (double *)du.data, *pv = (double *)dv.data, tmp;

    res.setTo(0);
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            io = r * istep + c * channels;
            fo = r * fstep + c;
            
            for (int k = 0; k < channels; k++)
            {
                if (LapPara[k] < 1E-20)
                    continue;
                
                tmp = pt[io+k] + px[io+k]*pu[fo] + py[io+k]*pv[fo];
                p[io+k] = psi_d(tmp * tmp);
            }
        }
    }
}

// calculate res = phi_1st(u^2 + v^2)
void OpticalFlow::phi_d(const Mat &u, const Mat &v, Mat &res)
{
    assert(matchAll(u, v) && matchAll(v, res) && u.type() == CV_64F &&
           res.type() == CV_64F);
    
    int rows = u.rows, cols = u.cols;
    int step = u.step / sizeof(double);
    double *p = (double *)res.data;

    Mat ux = dx(u), uy = dy(u), vx = dx(v), vy = dy(v);
    double *pux = (double *)ux.data, *puy = (double *)uy.data;
    double *pvx = (double *)vx.data, *pvy = (double *)vy.data;
    
    int r, c, o;
    double tmp;
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            o = r * step + c;
            tmp = pux[o]*pux[o] + puy[o]*puy[o] + pvx[o]*pvx[o] + pvy[o]*pvy[o];
            p[o] = phi_d(tmp);            
        }
    }
}

// get Ix, Iy, It
// smooth im1 -> sim1, im2 -> sim2
// im = 0.4 * sim1 + 0.6 * sim2
// Ix = dx(im), Iy = dy(im), It = sim2 - sim1
void OpticalFlow::getGrads(const Mat &im1, const Mat &im2, Mat &Ix, Mat &Iy, Mat &It)
{
    Mat sIm1, sIm2;
    CvSize size = cvSize(5, 5);
    
    GaussianBlur(im1, sIm1, size, 0);
    GaussianBlur(im2, sIm2, size, 0);

    Mat im = 0.4 * sIm1 + 0.6 * sIm2;
    Ix = dx(im);
    Iy = dy(im);
    It = sIm2 - sIm1;

    int rows = im1.rows, cols = im1.cols;
    Mat cix(rows, cols, CV_64F), ciy(rows, cols, CV_64F), cit(rows, cols, CV_64F);
}

void OpticalFlow::sanityCheck(const Mat &Ix, const Mat &Iy, const Mat &It,
                              const Mat &du, const Mat &dv)
{
    assert(matchAll(Ix, Iy) && matchAll(Iy, It) && matchAll(du, dv) &&
           Ix.depth() == CV_64F && du.type() == CV_64F);

    int rows = Ix.rows, cols = Ix.cols, channels = Ix.channels();
    int istep = Ix.step / sizeof(double);
    int fstep = du.step / sizeof(double);
    double *pix = (double *)Ix.data, *piy = (double *)Iy.data, *pit = (double *)It.data;
    double *pdu = (double *)du.data, *pdv = (double *)dv.data;
    int fo, io;
    double error = 0;

    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            fo = r * fstep + c;
            for (int k = 0; k < channels; k++)
            {
                io = r * istep + c * channels + k;
                error += fabs(pix[io]*pdu[fo] + piy[io]*pdv[fo] + pit[io]);
            }
        }
    }

    error = error / (rows * cols);
    printf("mean error for |Ix*u + Iy*v + It| = %.2f\n", error);
}

void OpticalFlow::estLapNoise(const Mat &im1, const Mat &im2)
{
    assert(matchAll(im1, im2));

    int rows = im1.rows, cols = im1.cols, channels = im1.channels(), offset;
    int step = im1.step / sizeof(double);
    double *p1 = (double *)im1.data, *p2 = (double *)im2.data;
    std::vector<int> total(channels, 0);
    double tmp;
    
    LapPara.resize(channels, 0);

    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            for (int k = 0; k < channels; k++)
            {
                tmp = fabs(p1[offset+k] - p2[offset+k]);
                if (tmp >= 0 && tmp < 1000000)
                {
                    LapPara[k] += tmp;
                    total[k]++;
                }
            }
        }
    }
    
    for (int k = 0; k < channels; k++)
    {
        if (total[k] == 0)
        {
            cout << "All the pixels are invalid in estimation Laplacian noise!!!" << endl;
			cout << "Something severely wrong happened!!!" << endl;
			LapPara[k] = 0.001;
		}
		else
			LapPara[k] /= total[k];
    }
}
