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
void OpticalFlow::warpImage(Mat &im1, Mat &im2, Mat &u, Mat &v, Mat &warp)
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

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            ioffset = r * istep + c * channels;
            foffset = r * fstep + c;
            nx = c + pu[foffset];
            ny = r + pv[foffset];
                
            if (nx < 0 || nx >= cols || ny < 0 || ny >= rows)
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
void OpticalFlow::compute(Mat &im1, Mat &im2, Mat &warp, Mat &u, Mat &v, double a_s,
                          int nInIter, int nOutIter, int nSORIter)
{
    int rows = im1.rows;
    int cols = im1.cols;
    int channels = im1.channels();

    Mat du(rows, cols, CV_64F), dv(rows, cols, CV_64F);
    Mat uu(rows, cols, CV_64F), vv(rows, cols, CV_64F);
    
    Mat It, Ix, Iy;
    Mat Ixy(rows, cols, CV_64FC(channels)), Ix2(rows, cols, CV_64FC(channels));
    Mat Iy2(rows, cols, CV_64FC(channels));
    Mat Ixt(rows, cols, CV_64FC(channels)), Iyt(rows, cols, CV_64FC(channels));
     
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
    int step = du.step / sizeof(double);
    int offset, tmp;

    // outer fixed point iteration for u, v
    for (int oiter = 0; oiter < nOutIter; oiter++)
    {
        cout << "outer fixed point iterations " << oiter << endl;

        // compute Ix, Iy and It
        getGrads(im1, warp, Ix, Iy, It);

        // inner fixed point iteration for du, dv
        du.setTo(0), dv.setTo(0);
        for (int iiter = 0; iiter < nInIter; iiter++)
        {
            cout << "inner fixed point iterations " << iiter << endl;

            uu = u + du;
            vv = v + dv;
            phi_d(uu, vv, phi_1st);
            psi_d(Ix, Iy, It, du, dv, psi_1st);

            // components for linear system
            multiply(Ix, Ix, Ix2);
            Ix2.mul(psi_1st);
            collapse(Ix2, ix2);
            
            multiply(Iy, Iy, Iy2);
            Iy2.mul(psi_1st);
            collapse(Iy2, iy2);
            
            multiply(Ix, Iy, Ixy);
            Ixy.mul(psi_1st);
            collapse(Ixy, ixy);
            
            multiply(Ix, It, Ixt);
            Ixt.mul(psi_1st);
            collapse(Ixt, ixt);
            
            multiply(Iy, It, Iyt);
            Iyt.mul(psi_1st);
            collapse(Iyt, iyt);

            weighted_lap(u, phi_1st, lapU);
            addWeighted(lapU, a_s, ixt, -1, 0, ixt);
            
            weighted_lap(v, phi_1st, lapV);
            addWeighted(lapV, a_s, iyt, -1, 0, iyt);

            // SOR iteration
            du.setTo(0), dv.setTo(0);
            
            const double omega = 1.8;
            double l, l_du, l_dv;

            for (int siter = 0; siter < nSORIter; siter++)
            {
                cout << "SOR iterations " << siter << endl;
                
                for (int r = 0; r < rows; r++)
                {
                    for (int c = 0; c < cols; c++)
                    {
                        offset = r * step + c;
                        l_du = 0, l_dv = 0, l = 0;
                        
                        if (r > 0)
                        {
                            tmp = offset - step;
                            l_du += pphi[offset] * pdu[tmp];
                            l_dv += pphi[offset] * pdv[tmp];
                            l += pphi[offset];
                        }
                        if (r < rows-1)
                        {
                            tmp = offset + step;
                            l_du += pphi[tmp] * pdu[tmp];
                            l_dv += pphi[tmp] * pdv[tmp];
                            l += pphi[tmp];
                        }
                        if (c > 0)
                        {
                            tmp = offset - 1;
                            l_du += pphi[offset] * pdu[tmp];
                            l_dv += pphi[offset] * pdv[tmp];
                            l += pphi[offset];
                        }
                        if (c < cols-1)
                        {
                            tmp = offset + 1;
                            l_du += pphi[tmp] * pdu[tmp];
                            l_dv += pphi[tmp] * pdv[tmp];
                            l += pphi[tmp];
                        }

                        l *= a_s;
                        l_du *= a_s;
                        l_dv *= a_s;

                        // du
                        l_du -= pixt[offset];
                        pdu[offset] = (1-omega) * pdu[offset] + omega / (pix2[offset] + l) * (l_du - pixy[offset] * pdv[offset]);
                        
                        // dv
                        l_dv -= piyt[offset];
                        pdv[offset] = (1-omega) * pdv[offset] + omega / (piy2[offset] + l) * (l_dv - pixy[offset] * pdu[offset]);
                    }
                }
            }
        }
        u += du;
        v += dv;
        warpImage(im1, im2, u, v, warp);
    }
}

void OpticalFlow::psi_d(Mat &Ix, Mat &Iy, Mat &It, Mat &du, Mat &dv, Mat &res)
{
    int rows = Ix.rows;
    int cols = Ix.cols;
    int channels = Ix.channels();
    int istep = Ix.step / sizeof(double);
    double *px = (double *)Ix.data, *py = (double *)Iy.data, *pz = (double *)It.data;
    double *p = (double *)res.data;
    
    int fstep = du.step / sizeof(double);
    double *pu = (double *)du.data, *pv = (double *)dv.data;
    
    int r, c, k, io, fo;
    double tmp;
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            io = r * istep + c * channels;
            fo = r * fstep + c;
            
            for (k = 0; k < channels; k++)
            {
                tmp = pz[io+k] + px[io+k]*pu[fo] + py[io+k]*pv[fo];
                p[io+k] = psi_d(tmp * tmp);
            }
        }
    }
}

void OpticalFlow::phi_d(Mat &u, Mat &v, Mat &res)
{
    int rows = u.rows;
    int cols = u.cols;
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

void OpticalFlow::getGrads(Mat &im1, Mat &im2, Mat &Ix, Mat &Iy, Mat &It)
{
    Mat sIm1, sIm2;
    CvSize size = cvSize(5, 5);
    
    GaussianBlur(im1, sIm1, size, 0);
    GaussianBlur(im2, sIm2, size, 0);
    
    Ix = dx(sIm2);
    Iy = dy(sIm2);
    It = sIm2 - sIm1;
}
