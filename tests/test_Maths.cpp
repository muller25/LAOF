#include "gtest/gtest.h"

#include "Image.h"
#include "Maths.h"
#include "ImageIO.h"

TEST(TestMaths, randFill)
{
    IImage im(3, 3);
    randFill(im, 0, 1);
    ASSERT_TRUE(im.ptr() != NULL);

    // imprint(im);
    
    for (int i = 0; i < im.nElements(); ++i)
        EXPECT_TRUE((im[i] <= 1) && (im[i] >= 0));

    DImage dm(5, 5);
    randFill(dm, -1., 0.);
    ASSERT_TRUE(dm.ptr() != NULL);

    // imprint(dm);
    
    for (int i = 0; i < dm.nElements(); ++i)
        EXPECT_TRUE((dm[i] <= 0) && (dm[i] >= -1));
}

TEST(TestMaths, rectSum)
{
    DImage im(10, 10, 10), dst, dstBF;
    randFill(im, -10., 10.);

    rectSum(dst, im);
    rectSumBF(dstBF, im);
    
    EXPECT_TRUE(equal(dst, dstBF));
}

TEST(TestMaths, cross)
{
    // 2 X 3, 3 X 2
    IImage im1(3, 2), im2(2, 3), res;
    randFill(im1, 1, 3);
    randFill(im2, 1, 5);
    cross(res, im1, im2);

    // printf("im1:\n");
    // imprint(im1);
    // printf("im2:\n");
    // imprint(im2);
    // printf("res:\n");
    // imprint(res);
    
    EXPECT_EQ(2, res.nWidth());
    EXPECT_EQ(2, res.nHeight());
}

TEST(TestImageProcess, split)
{
    int height = 100, width = 100, channels = 10, offset;
    DImage im(width, height, channels);
    std::vector<DImage> vec;

    for (int i = 0; i < 10; ++i)
    {
        randFill(im, -1000., 1000.);
        vec.clear();
        split(vec, im);

        for (int h = 0; h < height; ++h)
        {
            for (int w = 0; w < width; ++w)
            {
                offset = h * width + w;
                for (int k = 0; k < channels; ++k)
                    EXPECT_EQ(im[offset*channels+k], vec[k][offset]);
            }
        }
    }
}
