#ifndef _ImageIO_H
#define _ImageIO_H

#include <cv.h>
#include <highgui.h>

#include "Image.h"

template <class T>
void imread(Image<T> &im, const char * filename)
{
    cv::Mat img = cv::imread(filename);
    
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
