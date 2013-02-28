#include "gtest/gtest.h"

#include "OpticalFlow.h"
#include "Utils.h"

#include <cstdlib>
#include <ctime>

OpticalFlow of;

// psi(|Ix * u + Iy * v + It|^2)
TEST(TestOpticalFlow, TestPsiD)
{
    Mat ix(3, 3, CV_64FC3), iy(3, 3, CV_64FC3), it(3, 3, CV_64FC3);
    Mat du(3, 3, CV_64F), dv(3, 3, CV_64F);
    Mat exp(3, 3, CV_64FC3), res(3, 3, CV_64FC3);

    srand(time(NULL));
    double ixVal, iyVal, itVal, duVal, dvVal, tmp;
    double low = -1000, up = 1000;
    const double range = up - low;
    for (int i = 0; i < 10; i++)
    {
        ixVal = low + (double)rand()/(RAND_MAX/range);
        iyVal = low + (double)rand()/(RAND_MAX/range);
        itVal = low + (double)rand()/(RAND_MAX/range);
        duVal = low + (double)rand()/(RAND_MAX/range);
        dvVal = low + (double)rand()/(RAND_MAX/range);
        
        ix.setTo(ixVal);
        iy.setTo(iyVal);
        it.setTo(itVal);
        du.setTo(duVal);
        dv.setTo(dvVal);
        tmp = ixVal * duVal + iyVal * dvVal + itVal;
        tmp *= tmp;
        exp.setTo(of.psi_d(tmp));
        of.psi_d(ix, iy, it, du, dv, res);
        EXPECT_TRUE(matrix_match<double>(exp, res));
    }
}

// phi(|ux^2 + uy^2 + vx^2 + vy^2|)
TEST(TestOpticalFlow, TestPhiD)
{
    double arr[][3] = {{1., 2., 3.},
                       {2., 3., 4.},
                       {3., 4., 5.}};
    Mat u(3, 3, CV_64F, arr), v(3, 3, CV_64F, arr);
    Mat exp(3, 3, CV_64F), res(3, 3, CV_64F);
    Mat ux = dx(u), uy = dy(u), vx = dx(v), vy = dy(v);

    ASSERT_TRUE(matrix_match<double>(ux, vx));
    ASSERT_TRUE(matrix_match<double>(uy, vy));
        
    Mat ux2, uy2, vx2, vy2, g;
    multiply(ux, ux, ux2);
    multiply(uy, uy, uy2);
    multiply(vx, vx, vx2);
    multiply(vy, vy, vy2);
    g = ux2 + uy2 + vx2 + vy2;
    
    int rows = g.rows, cols = g.cols, offset;
    int step = g.step / sizeof(double);
    double *pg = (double *)g.data, *pe = (double *)exp.data;
    
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c;
            pe[offset] = of.phi_d(pg[offset]);
        }
    }
    
    of.phi_d(u, v, res);
    EXPECT_TRUE(matrix_match<double>(exp, res));
}

TEST(TestOpticalFlow, TestWarpImage)
{
    int rows = 3, cols = 3;
    Mat im1(rows, cols, CV_64F), im2(rows, cols, CV_64F), warp(rows, cols, CV_64F);
    Mat u = Mat::zeros(rows, cols, CV_64F), v = Mat::zeros(rows, cols, CV_64F);
    Mat exp(rows, cols, CV_64F), tmp(rows, cols, CV_64F);
    
    RNG rng(time(NULL));
    rng.fill(im1, RNG::UNIFORM, 0, 255.);
    rng.fill(im2, RNG::UNIFORM, 0, 255.);
    
    of.warpImage(im1, im2, u, v, warp);
    EXPECT_TRUE(matrix_match<double>(im2, warp));

    srand(time(NULL));
    double low = -1, high = 1;
    double range = high - low;
    double uVal, vVal;
    
    for (int i = 0; i < 10; i++)
    {
        uVal = low + rand() / (RAND_MAX/range);
        vVal = low + rand() / (RAND_MAX/range);

        cout << "uVal: " << uVal << ", vVal: " << vVal << endl;
        u.setTo(uVal), v.setTo(vVal);
        if (uVal >= 0)
        {
            for (int c = 0; c < cols-1; c++)
                tmp.col(c) = im2.col(c+1) * uVal + im2.col(c) * (1-uVal);
            im1.col(cols-1).copyTo(tmp.col(cols-1));
        } else {
            for (int c = 1; c < cols; c++)
                tmp.col(c) = im2.col(c) * (1 - fabs(uVal)) + im2.col(c-1) * fabs(uVal);
            im1.col(0).copyTo(tmp.col(0));
        }

        if (vVal >= 0)
        {
            for (int r = 0; r < rows-1; r++)
                exp.row(r) = tmp.row(r+1) * vVal + tmp.row(r) * (1-vVal);
        } else {
            for (int r = 1; r < rows; r++)
                exp.row(r) = tmp.row(r) * (1 - fabs(vVal)) + tmp.row(r-1) * fabs(vVal);
        }

        if (uVal >= 0)
            im1.col(cols-1).copyTo(exp.col(cols-1));
        else
            im1.col(0).copyTo(exp.col(0));

        if (vVal >= 0)
            im1.row(rows-1).copyTo(exp.row(rows-1));
        else
            im1.row(0).copyTo(exp.row(0));
        
        of.warpImage(im1, im2, u, v, warp);
        
        EXPECT_TRUE(matrix_match<double>(exp, warp));
    }
}
