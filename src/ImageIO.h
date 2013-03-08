#ifndef _ImageIO_H
#define _ImageIO_H

#include <cv.h>
#include <highgui.h>

#include "Image.h"

template <class T>
void imread(Image<T> &im, const char *filename)
{
    printf("reading image %s\n", filename);
    
    cv::Mat img = cv::imread(filename);
    assert(img.data != NULL);
    
    if (im.isFloat())
    {
        UCImage uc;
        uc.convertFrom(img);
        im2double(im, uc);
        return;
    }

    im.convertFrom(img);
}

template <class T>
bool imreadf(Image <T> &im, const char *filename)
{
    cv::FileStorage fs(filename, cv::FileStorage::READ);
    cv::Mat m;

    if (!fs.isOpened())
    {
        fs.release();
        return false;
    }
    
    fs["matrix"] >> m;
    im.convertFrom(m);
    fs.release();
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
    std::cout << "[";
    for (int h = 0; h < height; ++h)
    {
        if (h > 0)
            std::cout << " ";

        std::cout << "[";
        
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            for (int k = 0; k < channels; ++k)
            {
                std::cout << p[offset+k];
                if (k < channels-1)
                    std::cout << ", ";
            }
            if (w < width-1)
                std::cout << "; ";
        }

        std::cout << "]";
        if (h < height-1)
            std::cout << "\n";
    }
    std::cout << "]\n";
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
