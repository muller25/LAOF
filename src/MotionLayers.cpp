#include "MotionLayers.h"

void MotionLayers::im2feature(DImage &features, const DImage &im)
{
    int size = im.nSize(), channels = im.nChannels(), offset;

    features.create(channels, size);
    for (int i = 0; i < size; ++i)
    {
        offset = i * channels;
        for (int k = 0; k < channels; ++k)
            features[offset+k] = im[offset+k];
    }
}

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

int MotionLayers::flowCluster(DImage &centers, DImage &layers, const DImage &im,
                              const DImage &flow, const DImage &rflow)
{
    assert(flow.match3D(rflow));
    
    const int start = 2;
    const int end = 10;
    const double na = 4;
    
    DImage features;
    UCImage labels;
    int clusters;

    // flow2feature(features, flow, rflow);
    im2feature(features, im);
    clusters = kmeans2(centers, labels, features, start, end, na,
                       MotionLayers::mdist);

    // count number of elements in each cluster
    int *count = new int[clusters];
    memset(count, 0, sizeof(int) * clusters);
    for (int i = 0; i < labels.nElements(); ++i)
        count[labels[i]]++;

    for (int i = 0; i < clusters; ++i)
        printf("label %d: %d\n", i, count[i]);

    delete []count;
    
    // change label map to layer image
    layers.create(flow.nWidth(), flow.nHeight());
    for (int i = 0; i < labels.nSize(); i++)
        layers[i] = (double)labels[i];
    
    return clusters;
}
