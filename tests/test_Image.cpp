#include "gtest/gtest.h"
#include "Image.h"
#include "ImageIO.h"
#include "Maths.h"

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

TEST(TestImage, isSymmetric)
{
    int n = 3;
    DImage dm(n, n);

    randFill(dm, -100., 100.);

    EXPECT_FALSE(dm.isSymmetric());
    
    // make dm symmetric
    for (int h = 0; h < n; ++h)
        for (int w = h+1; w < n; ++w)
            dm[h*n+w] = dm[w*n+h];

    EXPECT_TRUE(dm.isSymmetric());
}


TEST(TestImage, eigen)
{
    int n = 3;
    DImage dm(n, n), values, vectors, eigenVector, tmp;

    randFill(dm, -100., 100.);

    // make dm symmetric
    for (int h = 0; h < n; ++h)
        for (int w = h+1; w < n; ++w)
            dm[h*n+w] = dm[w*n+h];

    dm.eigen(values, vectors);            
    
    eigenVector.create(1, n);
    for (int i = 0; i < n; ++i)
    {
        for (int w = 0; w < n; ++w)
            eigenVector[w] = vectors[i*n+w];

        cross(tmp, dm, eigenVector);
        multiply(eigenVector, values[i]);
        
        EXPECT_TRUE(equal(tmp, eigenVector));
    }
}
