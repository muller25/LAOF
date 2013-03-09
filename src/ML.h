#ifndef _ML_H
#define _ML_H

#include "Image.h"
#include "Maths.h"

#include <ctime>

template<class T>
void initClusterCenters(Image<T> &centers, const Image<T> &samples, const int clusters,
                        double (*distance)(T*, T*, int)=dist2)
{
    printf("using kmeans++ to init cluster centers...");
    
    int width = sample.nWidth(), height = samples.nHeight(), idx;
    centers.create(width, clusters);

    T *pc = centers.ptr(), *ps = samples.ptr();
    double *dist = new double[height];
    double sum = 0;
    
    // find the first cluster center uniformly randomly
    srand(time(NULL));
    idx = rand() % height;
    memcpy(pc, ps+idx*width, width*sizeof(T));

    for (int h = 0; h < height; ++h)
    {
        dist[h] = (*distance)(ps+h*width, pc, width);
        sum += dist[h];
    }
    
    // find remain cluster centers
    for (int c = 1; c < clusters; ++c)
    {
        
    }

    delete []dist;
    printf("done\n");
}

// each row contains one sample
template <class T>
double kmeans(Image<T> &centers, UCImage &labels,
              const Image<T> &samples, const int clusters,
              double (*distance)(T *, T *, int)=dist2)
{
    printf("running kmeans...\n");
    assert(samples.nChannels() == 1 && clusters > 0);
    
    int width = samples.nWidth(), height = samples.nHeight(), idx;    
    const int maxIter = 1000;
    const double threshold = 0.001 * width;
    double best_compactness = DBL_MAX, max_shift = DBL_MAX, compactness, shift, dist, tmp;
    int *counters = new int[clusters];
    Image<T> nCenters(width, clusters);
    
    if (centers.isEmpty())
        initClusterCenters(centers, samples, clusters, distance);// kmeans++
    labels.create(1, height);
    T *ps = samples.ptr(), *pc = centers.ptr(), *pnc = nCenters.ptr();

    for (int iter = 0; iter < maxIter && max_shift > threshold; ++iter)
    {
        // calc data point to nearest cluster center
        printf("calc data point to nearest cluster center\n");
        nCenters.setTo(0);
        memset(counters, 0, sizeof(int) * clusters);
        for (int h = 0; h < height; ++h)
        {
            dist = DBL_MAX;
            idx = -1;
            for (int c = 0; c < clusters; ++c)
            {
                tmp = (*distance)(ps+h*width, pc+c*width, width);
                if (tmp < dist)
                {
                    dist = tmp;
                    idx = c;
                }
            }

            // assign label
            labels[h] = idx;
            counters[idx]++;
            for (int w = 0; w < width; ++w)
                pnc[idx*width+w] += ps[h*width+w];
        }

        // check whether there is cluster without element
        for (int c = 0; c < clusters; ++c)
        {
            if (counters[c] > 0) continue;

            printf("find cluster without element, create new cluster...");
            
            printf("done\n");
        }
        
        // calc new cluster centers && center shift
        max_shift = 0;
        for (int c = 0; c < clusters; ++c)
        {
            shift = 0;
            for (int w = 0; w < width; ++w)
            {
                idx = c * width + w;
                pnc[idx] = (double)pnc[idx] / counters[c];
                shift += fabs(pnc[idx] - pc[idx]);
            }

            if (max_shift < shift) max_shift = shift;
        }

        // calc compactness of clusters
        compactness = 0;
        for (int h = 0; h < height; ++h)
            compactness += dist2(ps+h*width, pc+labels[h]*width, width);

        if (compactness < best_compactness)
        {
            best_compactness = compactness;
            nCenters.copyTo(centers);
        }

        printf("iterations %d, compactness: %.6f, max_shift: %.6f\n", iter, best_compactness, max_shift);
    }

    delete []counters;
    return best_compactness;
}

#endif
