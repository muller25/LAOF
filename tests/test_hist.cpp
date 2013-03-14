#include "gtest/gtest.h"
#include "Image.h"
#include "Maths.h"
#include "Histogram.h"
#include "ImageIO.h"

TEST(TestHist, calcHist)
{
    int width = 4, height = 4, channels = 1;
    int minV = 0, maxV = 10;
    
    IImage m(width, height, channels), hist;
    Image<bool> mask(width, height, 1, 1);

    randFill(m, minV, maxV);
    calcHist(hist, m, mask, minV, maxV, 5);

    int total = 0;
    for (int i = 0; i < hist.nElements(); ++i)
        total += hist[i];
    
    EXPECT_EQ(mask.nonZeros(), total);

    mask.set(false);
    calcHist(hist, m, mask, minV, maxV, 5);

    total = 0;
    for (int i = 0; i < hist.nElements(); ++i)
        total += hist[i];

    EXPECT_EQ(0, total);
}

TEST(TestOMHist, calcOMHist)
{
    int width = 4, height = 4, channels = 2;
    double minV = -5, maxV = 5;
    
    DImage m(width, height, channels), hist;
    Image<bool> mask(width, height, 1, 1);

    randFill(m, minV, maxV);
    calcOMHist(hist, m, mask, 18);

    double total = 0;
    for (int i = 0; i < hist.nElements(); ++i)
        total += hist[i];

    double expect = 0, x, y;
    for (int i = 0; i < m.nSize(); ++i)
    {
        if (mask.isZero(i)) continue;

        x = m[i*2];
        y = m[i*2+1];
        expect += sqrt(x*x + y*y) * (atan2(y, x) + PI) / PI * 180.;
    }

    EXPECT_DOUBLE_EQ(expect, total);
}
