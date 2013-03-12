#include "gtest/gtest.h"
#include "Image.h"

TEST(TestImage, isEmpty)
{
    DImage dm;
    IImage im;
    
    EXPECT_TRUE(dm.isEmpty());

    dm.create(1, 1);
    EXPECT_FALSE(dm.isEmpty());

    dm.create(3, 3);
    EXPECT_FALSE(dm.isEmpty());
    
    dm.release();
    EXPECT_TRUE(dm.isEmpty());

    im.create(3, 3);
    im.copyTo(dm);
    EXPECT_FALSE(dm.isEmpty());

    im.release();
    EXPECT_TRUE(im.isEmpty());

    im.copyFrom(dm);
    EXPECT_FALSE(dm.isEmpty());
}
