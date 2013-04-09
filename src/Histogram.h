#ifndef _Histogram_H
#define _Histogram_H

#include "Image.h"

// normalized lab histogram
template<class I, class M>
void LabHist(DImage &hist, Image<T> &im, Image<M> &mask)
{
    assert(im.match2D(mask) && im.nChannels() == 3);

    // convert to lab color space
    cv::Mat imat, lab, maskMat, tmp;

    if (typeid(T) != typeid(float))
    {
        im.covertTo(tmp);
        tmp.convertTo(imat, CV_32FC3);
    }
    else
        im.convertTo(imat);
        
    cv::cvtColor(imat, lab, CV_BGR2Lab);

    // convert mask
    if (typeid(M) != typeid(uchar))
    {
        mask.convertTo(tmp);
        tmp.convertTo(maskMat, CV_8U);
    }
    else
        mask.convertTo(maskMat);

    const int histSize = {10, 12, 12};
    const int channels = {0, 1, 2};
                          
    // l varies from [0, 100]
    const float lranges[] = {0, 101};
    // a varies from [-127, 127]
    const float aranges[] = {-127, 128};
    // b varies from [-127, 127]
    const float branges[] = {-127, 128};

    const float *ranges[] = {lranges, aranges, branges};

    cv::Mat histMat;
    
                            
}

#endif
