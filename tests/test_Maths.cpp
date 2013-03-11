#include "gtest/gtest.h"

#include "Image.h"
#include "Maths.h"
#include "ImageIO.h"

TEST(TestRandFill, randFill)
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

TEST(TestRectSum, rectSum)
{
    IImage im(5, 5), dst, dstBF;
    randFill(im, 0, 2);

    rectSum(dst, im);
    rectSumBF(dstBF, im);

    imprint(im);
    printf("O(1)\n");
    imprint(dst);

    printf("brute force\n");
    imprint(dstBF);
    
    EXPECT_TRUE(equal(dst, dstBF));
}
