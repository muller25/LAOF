#include "gtest/gtest.h"

#include "Maths.h"
#include "Utils.h"

template<class T, int channels>
void tests()
{
    int depth = DataType<T>::depth;
    Mat mat(2, 2, CV_MAKETYPE(depth, channels));
    T data1[4] = {1, 2, 3, 4};
    T data2[8] = {1, 1, 2, 2, 3, 3, 4, 4};
    T data3[12] = {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4};
    T res[channels], exp[channels];
    
    switch (channels)
    {
    case 1:
        setData(mat, data1);
        break;
    case 2:
        setData(mat, data2);
        break;
    case 3:
        setData(mat, data3);
        break;
    default:
        assert(channels <= 3);
        break;
    }

    //   std::cout << mat << std::endl;
    
    std::fill(exp, exp+channels, 2.5);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 0.5, 0.5, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 1);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 0, 0, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 2);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 1, 0, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 3);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 0, 1, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 4);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 1, 1, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 1.5);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 0.5, 0, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 2.8);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 0, 0.9, res);
    EXPECT_TRUE(array_match(exp, res, channels));

    std::fill(exp, exp+channels, 2.65);
    std::fill(res, res+channels, 0);
    biInterpolate(mat, 0.33, 0.66, res);
    EXPECT_TRUE(array_match(exp, res, channels));
}

TEST(TestBiInterpolation, MultiChannel)
{
    tests<uchar, 1>();
    tests<uchar, 2>();
    tests<uchar, 3>();

    tests<double, 1>();
    tests<double, 2>();
    tests<double, 3>();
}

TEST(TestGrads, MultiChannel)
{
    // single channel
    double arr[][3] = {{1., 2., 3.},
                       {2., 3., 4.},
                       {3., 4., 5.}};
    double exp_1st[][3] = {{0.5, 14./12, 0.5},
                           {0.5, 14./12, 0.5},
                           {0.5, 14./12, 0.5}};
    Mat m(3, 3, CV_64F, arr);
    Mat exp_dx(3, 3, CV_64F, exp_1st);

    Mat gx, gy;

    gradX(m, gx);
    EXPECT_TRUE(matrix_match<double>(exp_dx, gx));

    gradY(m, gy);
    EXPECT_TRUE(matrix_match<double>(exp_dx.t(), gy));

    double exp_2st[][3]={{14./12, 0, -14./12},
                         {14./12, 0, -14./12},
                         {14./12, 0, -14./12}};
    Mat exp_dxx(3, 3, CV_64F, exp_2st);

    Mat gxx = gradXX(m);
    EXPECT_TRUE(matrix_match<double>(exp_dxx, gxx));

    Mat gyy = gradYY(m);
    EXPECT_TRUE(matrix_match<double>(exp_dxx.t(), gyy));

    double exp_2st1[][3] = {{0, 0, 0},
                            {0, 0, 0},
                            {0, 0, 0}};
    Mat exp_dxy(3, 3, CV_64F, exp_2st1);

//    Mat gxy = gradXY(m);
//    EXPECT_TRUE(matrix_match<double>(exp_dxy, gxy));

    // multi-channel
    double arr3[][9] = {{1., 1., 1., 2., 2., 2., 3., 3., 3.},
                        {2., 2., 2., 3., 3., 3., 4., 4., 4.},
                        {3., 3., 3., 4., 4., 4., 5., 5., 5.}};
    double exp3_1st[][9] = {{0.5, 0.5, 0.5, 14./12, 14./12, 14./12, 0.5, 0.5, 0.5},
                            {0.5, 0.5, 0.5, 14./12, 14./12, 14./12, 0.5, 0.5, 0.5},
                            {0.5, 0.5, 0.5, 14./12, 14./12, 14./12, 0.5, 0.5, 0.5}};
    Mat m3(3, 3, CV_64FC3, arr3);
    Mat exp3_dx(3, 3, CV_64FC3, exp3_1st);

    gradX(m3, gx);
    EXPECT_TRUE(matrix_match<double>(exp3_dx, gx));

    gradY(m3, gy);
    EXPECT_TRUE(matrix_match<double>(exp3_dx.t(), gy));

    double exp3_2st[][9]={{14./12, 14./12, 14./12, 0, 0, 0, -14./12, -14./12, -14./12},
                          {14./12, 14./12, 14./12, 0, 0, 0, -14./12, -14./12, -14./12},
                          {14./12, 14./12, 14./12, 0, 0, 0, -14./12, -14./12, -14./12}};
    Mat exp3_dxx(3, 3, CV_64FC3, exp3_2st);

    gxx = gradXX(m3);
    EXPECT_TRUE(matrix_match<double>(exp3_dxx, gxx));

    gyy = gradYY(m3);
    EXPECT_TRUE(matrix_match<double>(exp3_dxx.t(), gyy));

    double exp3_2st1[][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    Mat exp3_dxy(3, 3, CV_64FC3, exp3_2st1);
    
//    gxy = gradXY(m3);
//    EXPECT_TRUE(matrix_match<double>(exp3_dxy, gxy));
}

TEST(TestCollapse, MultiChannel)
{
    int arr3[][9] = {{1, 2, 3, 4, 5, 6, 7, 8, 9},
                     {2, 3, 4, 5, 6, 7, 8, 9, 10},
                     {3, 4, 5, 6, 7, 8, 9, 10, 11}};
    int arr[][3] = {{2, 5, 8},
                    {3, 6, 9},
                    {4, 7, 10}};
    Mat m3(3, 3, CV_32SC3, arr3), m1(3, 3, CV_32S, arr);
    Mat exp(3, 3, CV_32S, arr), res(3, 3, CV_32S);

    collapse<int>(m1, res);
    EXPECT_TRUE(matrix_match<int>(exp, res));

    collapse<int>(m3, res);
    EXPECT_TRUE(matrix_match<int>(exp, res));

    double darr3[][9] = {{.1, .2, .3, .4, .5, .6, .7, .8, .9},
                         {.2, .3, .4, .5, .6, .7, .8, .9, .1},
                         {.3, .4, .5, .6, .7, .8, .9, .1, .11}};
    double darr[][3] = {{.2, .5, .8},
                        {.3, .6, .6},
                        {.4, .7, .37}};
    Mat dm1(3, 3, CV_64F, darr), dm3(3, 3, CV_64FC3, darr3);
    Mat dexp(3, 3, CV_64F, darr), dres(3, 3, CV_64F);

    collapse<double>(dm1, dres);
    EXPECT_TRUE(matrix_match<double>(dexp, dres));

    collapse<double>(dm3, dres);
    EXPECT_TRUE(matrix_match<double>(dexp, dres));
}

TEST(TestWeightedLap, Lap2D)
{
    double arr[][3] = {{1., 2., 3.},
                       {2., 3., 4.},
                       {3., 4., 5.}};
    double earr[][3] = {{-2., -3., -1.},
                        {-3., -2., 2.},
                        {-1., 2., 8.}};
    Mat flow(3, 3, CV_64F, arr), weight(3, 3, CV_64F, arr), dst(3, 3, CV_64F);
    Mat exp(3, 3, CV_64F, earr);
            
    weighted_lap(flow, weight, dst);
    EXPECT_TRUE(matrix_match<double>(exp, dst));
}

TEST(TestWeightedLap, Lap3D)
{
}
