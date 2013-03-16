#ifndef _MotionLayers_H
#define _MotionLayers_H

#include "Image.h"

class MotionLayers
{
public:
    void spatialInfo(DImage &info, const DImage &im);
    void imInfo(DImage &info, const DImage &im);
    void flowInfo(DImage &info, const DImage &flow);

    int cluster(DImage &centers, DImage &layers, const DImage &im,
                const DImage &flow, const DImage &rflow,
                int start=2, int end=10, double na=15, bool reArrange=false);

    int cluster(DImage &centers, DImage &layers,
                const DImage &features, int width, int height,
                int start=2, int end=10, double na=15, bool reArrange=false);

    void refine(DImage &centers, DImage &layers, int labels,
                const DImage &im1, const DImage &im2,
                const DImage &features);

    // make sure layer order
    template <class T>
    void reArrangeLabels(Image<T> &layers, int labels);

    // create cluster centers by labels
    template <class T, class T1>
    void createCenterByLabels(Image<T> &centers, int numOfLabels,
                              const Image<T1> &labels, const Image<T> &samples);
    
   
    static double smoothFn(int p1, int p2, int l1, int l2, void *pData);
    static double mydist(double *p1, double *p2, int start, int end);

private:
    const static int sWidth = 2;// spatial info width
    const static int iWidth = 3;// image info width
    const static int fWidth = 2;// flow info width
};

// arrange layers from largest component to smallest
template <class T>
void MotionLayers::reArrangeLabels(Image<T> &layers, int labels)
{
    assert(layers.nChannels() == 1);
    
    int *count = new int[labels];
    int *mapping = new int[labels];
    int idx, maxNum, maxid, id;

    memset(count, 0, sizeof(int) * labels);
    for (int i = 0; i < layers.nElements(); ++i)
    {
        idx = layers[i] + 0.5;
        count[idx]++;
    }

    id = 0, maxid = -1;
    for (int i = 0; i < labels; ++i)
    {
        maxNum = -1;
        for (int j = 0; j < labels; ++j)
        {
            if (maxNum < count[j]) {
                maxNum = count[j];
                maxid = j;
            }
        }

        mapping[maxid] = id++;
        count[maxid] = -1;
    }

    for (int i = 0; i < layers.nElements(); ++i)
    {
        idx = layers[i] + 0.5;
        layers[i] = mapping[idx];
    }
    
    delete []mapping;
    delete []count;
}

template <class T, class T1>
void MotionLayers::createCenterByLabels(Image<T> &centers, int numOfLabels,
                                        const Image<T1> &labels, const Image<T> &samples)
{
    // printf("create center by label... ");
    
    int *count = new int[numOfLabels];
    int width = samples.nWidth(), lsize = labels.nSize(), idx;
    
    centers.create(width, numOfLabels);
    memset(count, 0, sizeof(int)*numOfLabels);

    for (int i = 0; i < lsize; ++i)
    {
        idx = labels[i];
        count[idx]++;
        for (int w = 0; w < width; ++w)
            centers[idx*width+w] += samples[i*width+w];
    }

    for (int i = 0; i < numOfLabels; ++i)
        for (int w = 0; w < width; ++w)
            centers[i*width+w] /= (double)count[i];
    
    delete []count;
    // printf("done\n");
}


#endif
