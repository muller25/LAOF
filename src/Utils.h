#ifndef _Utils_H
#define _Utils_H

#include "Image.h"

// 根据分层id生成分层
template <class L, class M>
void genLayerMask(Image<M> &mask, const Image<L> &layers, int layerID)
{
    mask.create(layers.nWidth(), layers.nHeight());
    for (int i = 0; i < layers.nElements(); ++i)
        if (fabs(layers[i] - layerID) < 0.5)
            mask[i] = 1;
}

// 用灰度图显示Layers
template<class L>
void showLayers(UCImage &img, const Image<L> &layers)
{
    assert(layers.nChannels() == 1);
    
    int maxVal = layers.max();
    int scale = 255 / (maxVal + 1);
    
    img.create(layers.nWidth(), layers.nHeight());
    for (int i = 0; i < layers.nElements(); ++i)
        img[i] = (int)(layers[i] + 1) * scale;
}

#endif
