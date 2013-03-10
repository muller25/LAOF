#include "MotionLayers.h"

// extract spatial info from image
void MotionLayers::spatialInfo(DImage &info, const DImage &im)
{
    int width = im.nWidth(), height = im.nHeight(), offset;
    info.create(sWidth, width*height);

    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * sWidth;
            info[offset] = (double)h / (height-1);
            info[offset+1] = (double)w / (width-1);
        }
    }
}

// extract intensity info from image
void MotionLayers::imInfo(DImage &info, const DImage &im)
{
    int size = im.nSize();

    DImage gray;
    desuarate(gray, im);
    
    info.create(iWidth, size);
    for (int i = 0; i < size; ++i)
        info[i] = gray[i];
}

// extract motion strength from flow
void MotionLayers::flowInfo(DImage &info, const DImage &flow)
{
    assert(flow.nChannels() == fWidth);
    int size = flow.nSize(), channels = flow.nChannels(), offset;
    DImage nf;
    
    normalizeChannels(nf, flow); // normalize is better
    // flow.copyTo(nf);

    info.create(channels, size);
    for (int i = 0; i < size; ++i)
    {
        offset = i * channels;
        for (int k = 0; k < channels; k++)
            info[offset+k] = nf[offset+k];
    }
}

int MotionLayers::cluster(DImage &centers, DImage &layers, const DImage &im,
                          const DImage &flow, const DImage &rflow)
{
    const int start = 2;
    const int end = 10;
    const double na = 10;
    
    DImage features;
    UCImage labels;
    int clusters;
    DImage fInfo, rfInfo;
    std::vector<DImage> vec;

    // DImage sInfo;
    // spatialInfo(sInfo, im);
    // vec.push_back(sInfo);

    // DImage iInfo;
    // imInfo(iInfo, im);
    // vec.push_back(iInfo);

    flowInfo(fInfo, flow);
    vec.push_back(fInfo);

    flowInfo(rfInfo, rflow);
    vec.push_back(rfInfo);
    
    mergew(features, vec);
    clusters = kmeans2(centers, labels, features, start, end, na,
                       MotionLayers::mydist);

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
