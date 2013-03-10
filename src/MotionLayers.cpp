#include "MotionLayers.h"

void MotionLayers::flow2feature(DImage &features,
                                const DImage &flow, const DImage &rflow)
{
    assert(flow.match3D(rflow) && flow.nChannels() == 2);
    
    int size = flow.nSize(), offset, foffset;
    DImage nf, nrf;
    
    normalizeChannels(nf, flow);
    normalizeChannels(nrf, rflow);
//    flow.copyTo(nf);
//    rflow.copyTo(nrf);
    features.create(4, size);
    for (int i = 0; i < size; ++i)
    {
        offset = i * 2;
        foffset = i * 4;

        features[foffset] = nf[offset];
        features[foffset+1] = nf[offset+1];
        features[foffset+2] = nrf[offset];
        features[foffset+3] = nrf[offset+1];
    }
}

int MotionLayers::flowCluster(DImage &centers, DImage &layers,
                              const DImage &flow, const DImage &rflow)
{
    assert(flow.match3D(rflow));
    
    const int start = 2;
    const int end = 20;
    const double na = 4;
    
    DImage features;
    UCImage labels;
    int clusters;
    
    flow2feature(features, flow, rflow);
    clusters = kmeans2(centers, labels, features, start, end, na,
                       MotionLayers::mdist);

    // change label map to layer image
    layers.create(flow.nWidth(), flow.nHeight());
    for (int i = 0; i < labels.nSize(); i++)
        layers[i] = (double)labels[i];
    
    return clusters;
}
