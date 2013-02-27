#include "gtest/gtest.h"
#include "Maths.h"

#include <cv.h>
using namespace cv;

#include <iostream>

template<class T>
bool array_match(T *exp, T *act, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (fabs(exp[i] - act[i]) <= 10E-6) continue;

        std::cout << "expected[" << i << "(" << exp[i]
             << ") != actual[" << i << "(" << act[i] << ")\n";
        return false;
    }

    return true;
}

template<class T>
void setData(Mat &m, const T *data)
{
    int rows = m.rows;
    int cols = m.cols;
    int channels = m.channels();
    int r, c, k, offset;
    int step = m.step / sizeof(T);
    T *ptr = (T *)m.data;
    
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            for (k = 0; k < channels; k++)
                ptr[offset+k] = data[offset+k];
        }
    }
}

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
