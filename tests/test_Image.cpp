#include "gtest/gtest.h"
#include "Image.h"
#include "ImageIO.h"
#include "Maths.h"

#include <cv.h>

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

TEST(TestImage, convertTo)
{
    cv::Mat mat;
    int width = 300, height = 300, offset, moffset, step;
    DImage dm(width, height);

    randFill(dm, -100., 100.);
    dm.convertTo(mat);

    step = mat.step1();
    double *dptr = (double *)mat.data;
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            moffset = h * step + w;
            EXPECT_DOUBLE_EQ(dm[offset], dptr[moffset]);
        }
    }

    Image<float> fm(width, height);

    randFill(fm, (float)-100., (float)100.);
    fm.convertTo(mat);

    step = mat.step1();
    float *fptr = (float *)mat.data;
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            moffset = h * step + w;
            EXPECT_DOUBLE_EQ(fm[offset], fptr[moffset]);
        }
    }
}

TEST(TestImage, convertFrom)
{
    int width = 300, height = 300, channels = 3, offset, moffset, step;
    cv::Mat mat(height, width, CV_8UC(channels));
    cv::RNG rng;
    DImage dm;
    
    rng.fill(mat, cv::RNG::UNIFORM, -100, 100);
    step = mat.step1();
    uchar *ptr = mat.data;

    dm.convertFrom(mat);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            moffset = h * step + w * channels;
            for (int k = 0; k < channels; ++k)
            {
                double tmp = ptr[moffset+k];
                EXPECT_DOUBLE_EQ(dm[offset+k], tmp);
            }
        }
    }

    IImage im;
    im.convertFrom(mat);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            moffset = h * step + w * channels;

            for (int k = 0; k < channels; ++k){
                int tmp = ptr[moffset+k];
                EXPECT_EQ(im[offset+k], tmp);
            }
        }
    }

    cv::Mat dmat(height, width, CV_64FC(channels));
    rng.fill(dmat, cv::RNG::UNIFORM, -100., 100.);
    step = dmat.step1();
    double *dptr = (double *)dmat.data;

    im.convertFrom(dmat);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            moffset = h * step + w * channels;
            for (int k = 0; k < channels; ++k)
            {
                int tmp;
                if (dptr[moffset + k] < 0)
                    tmp = dptr[moffset + k] - 0.5;
                else
                    tmp = dptr[moffset + k] + 0.5;
                
                EXPECT_EQ(im[offset+k], tmp);
            }
        }
    }
}
