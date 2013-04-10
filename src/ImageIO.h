#ifndef _ImageIO_H
#define _ImageIO_H

#include <cv.h>
#include <highgui.h>

#include "Image.h"
#include "ImageProcess.h"

enum ImageType{GRAY, COLOR};
    
template <class T>
void imread(Image<T> &im, const char *filename, ImageType it=COLOR)
{
    printf("reading image %s... ", filename);
    cv::Mat img;
    
    switch (it)
    {
    case GRAY:
        img = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
        break;
    case COLOR:
        img = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
        break;
    }
    
    assert(img.data != NULL);
    
    if (im.isFloat())
    {
        UCImage uc;
        uc.convertFrom(img);
        im2double(im, uc);
        printf("done\n");
        return;
    }

    im.convertFrom(img);
    printf("done\n");
}

template <class T>
bool imreadf(Image <T> &im, const char *filename)
{
    printf("reading file %s... ", filename);
    
    cv::FileStorage fs(filename, cv::FileStorage::READ);
    cv::Mat m;

    if (!fs.isOpened())
    {
        fs.release();
        printf("failed\n");
        return false;
    }
    
    fs["matrix"] >> m;
    im.convertFrom(m);
    fs.release();
    printf("done\n");
    return true;
}

template <class T>
void imwrite(const char *filename, const Image<T> &im)
{
    cv::Mat img;

    if (im.isFloat())
    {
        cv::Mat tmp;
        im.convertTo(tmp);
        tmp.convertTo(img, CV_8U, 255);
    }
    else im.convertTo(img);
    
    cv::imwrite(filename, img);
}

template <class T>
void imwritef(const char *filename, const Image<T> &im)
{
    cv::FileStorage fs(filename, cv::FileStorage::WRITE);
    cv::Mat m;
    im.convertTo(m);
    fs << "matrix" << m;
    fs.release();
}

template <class T>
void imprint(const Image<T> &im)
{
    T *p = im.ptr();
    int width = im.nWidth(), height = im.nHeight(), channels = im.nChannels();
    int offset;
    bool isFloat = im.isFloat();
    const char *iFmt = "%3d";
    const char *dFmt = "%.2f";
    
    printf("[");
    for (int h = 0; h < height; ++h)
    {
        if (h > 0)
            printf(" ");

        printf("[");
        
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            for (int k = 0; k < channels; ++k)
            {
                if (isFloat) printf(dFmt, p[offset+k]);
                else printf(iFmt, p[offset+k]);
                if (k < channels-1)
                    printf(", ");
            }
            if (w < width-1)
                printf("; ");
        }

        printf("]");
        if (h < height-1)
            printf("\n");
    }
    printf("]\n");
}

template <class T>
void imshow(const char *title, const Image<T> &im)
{
    cv::Mat img;

    im.convertTo(img);
    cv::imshow(title, img);
}

template <class T>
void imwait(T code)
{
    cv::waitKey(code);
}

#endif
