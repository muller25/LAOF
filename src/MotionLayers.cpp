#include "GCoptimization.h"
#include "MotionLayers.h"
#include "ML.h"
#include "Maths.h"
#include "ImageProcess.h"

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
            info[offset] = h;
            info[offset+1] = w;
        }
    }
}

// extract intensity info from image
void MotionLayers::imInfo(DImage &info, const DImage &im)
{
    int size = im.nSize(), channels = im.nChannels();
    assert(iWidth == channels);
    
    info.create(channels, size);
    for (int i = 0; i < im.nElements(); ++i)
        info[i] = im[i];
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
                          const DImage &flow, const DImage &rflow,
                          int start, int end, double na, bool reArrange)
{
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

    // rearrange labels
    if (reArrange)
    {
        printf("re-arrange labels...\n");
        reArrangeLabels(labels, clusters);
        createCenterByLabels(centers, clusters, labels, features);
    }

    // change label map to layer image
    layers.create(flow.nWidth(), flow.nHeight());
    for (int i = 0; i < labels.nSize(); i++)
        layers[i] = (double)labels[i];
    
    return clusters;
}

int MotionLayers::cluster(DImage &centers, DImage &layers, const DImage &features,
                          int width, int height, int start, int end, double na, bool reArrange)
{
    printf("running kmeans...\n");
    int clusters;
    UCImage labels;
    
    clusters = kmeans2(centers, labels, features, start, end, na,
                       MotionLayers::mydist);

    // rearrange labels
    if (reArrange)
    {
        printf("re-arrange labels...\n");
        reArrangeLabels(labels, clusters);
        createCenterByLabels(centers, clusters, labels, features);
    }
    
    // change label map to layer image
    layers.create(width, height);
    for (int i = 0; i < labels.nSize(); i++)
        layers[i] = (double)labels[i];

    return clusters;
}

void MotionLayers::refine(DImage &centers, DImage &layers, int labels,
                          const DImage &im1, const DImage &im2,
                          const DImage &features)
{
    assert(centers.ptr() != NULL);
    
    int size = im1.nSize(), cwidth = centers.nWidth();
    int width = im1.nWidth(), height = im1.nHeight(), channels = im1.nChannels();
    double *data = new double[labels*size];
    DImage extra;
    std::vector<DImage> vec;

    vec.push_back(im1);
    vec.push_back(im2);
    mergec(extra, vec);

    // set data term
    for (int i = 0; i < size; ++i)
    {
        for (int l = 0; l < labels; ++l)
        {
            data[i*labels+l] = mydist(features.ptr()+i*cwidth, centers.ptr()+l*cwidth, 0, cwidth);
            data[i*labels+l] += dist2(im1.ptr()+i*channels, im2.ptr()+i*channels, 0, channels);
        }
    }
    
    try{
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width, height, labels);
        gc->setDataCost(data);
        gc->setSmoothCost(MotionLayers::smoothFn, extra.ptr());
        
        printf("Before optimization energy is %.6f\n",gc->compute_energy());
        // gc->expansion(2);
        gc->swap(2);
        printf("After optimization energy is %.6f\n",gc->compute_energy());

        for (int  i = 0; i < size; i++)
            layers[i] = gc->whatLabel(i);

        // rearrange labels
        // reArrangeLabels(layers, labels);

        createCenterByLabels(centers, labels, layers, features);
        
		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

    delete []data;
}

double MotionLayers::smoothFn(int p1, int p2, int l1, int l2, void *pData)
{
    const int dataWidth = 6;
    const double weight = 1;
    const double penalty = 1;
    const double sigma = 2.5;
    double *ptr = (double *)pData;
    double cost;
    
    if (l1 == l2) return 0;

    cost = dist1(ptr+p1*dataWidth, ptr+p2*dataWidth, 0, dataWidth);
    cost = weight * exp(-cost*cost/(2*sigma*sigma)) + penalty;

    return cost;
}

double MotionLayers::mydist(double *p1, double *p2, int start, int end)
{
    double dist = 0;
    int idx = 0;
    
    // spatial info
    // dist += dist2(p1, p2, idx, idx+sWidth);
    // idx += sWidth;
    
    // image info
    // dist += dist2(p1, p2, idx, idx+iWidth);
    // idx += iWidth;
    
    // flow info
    dist += dist2(p1, p2, idx, idx+fWidth);
    dist += (1 - (similarity(p1, p2, idx, idx+fWidth) + 1) / 2);
    idx += fWidth;
    
    // reverse flow info
    dist += dist2(p1, p2, idx, idx+fWidth);
    dist += (1 - (similarity(p1, p2, idx, idx+fWidth) + 1) / 2);
    
    return dist;
}
