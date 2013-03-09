#ifndef _ML_H
#define _ML_H

#include "Image.h"
#include "Maths.h"

#include <ctime>

template<class T>
void findClusterCenters(Image<T> &centers, const Image<T> &samples, int clusters,
                        double *dist, double sum, int *target, int need2deter,
                        double (*distance)(T*, T*, int))
{
    // printf("find cluster centers...");
    
    const int trials = 3;
    int width = samples.nWidth(), height = samples.nHeight();
    int curIdx, bestIdx, h;
    T *ps = samples.ptr(), *pc = centers.ptr();
    double curSum, bestSum, tmp;
    double *curDist = new double[height];
    double *bestDist = new double[height];
        
    srand(time(NULL));
    for (int i = 0; i < need2deter; ++i)
    {
        // try to find cluster centers
        bestSum = DBL_MAX, bestIdx = -1;
        for (int t = 0; t < trials; ++t)
        {
            // proportional random
            curSum = sum * (double)rand() / RAND_MAX;
            for (curIdx = 0; curIdx < height; ++curIdx)
                if ((curSum -= dist[curIdx]) <= 0) break;

            assert(curIdx < height);
            
            // find minimum distance to cluster centers
            for (curSum = 0, h = 0; h < height; ++h)
            {
                tmp = (*distance)(ps+h*width, ps+curIdx*width, width);
                curDist[h] = std::min(dist[h], tmp);
                curSum += curDist[h];
            }

            // is this one better than previous one?
            if (curSum < bestSum)
            {
                bestSum = curSum;
                bestIdx = curIdx;
                memcpy(bestDist, curDist, sizeof(double)*height);
            }
        }
        
        assert(bestIdx >= 0);
        sum = bestSum;
        memcpy(dist, bestDist, sizeof(double)*height);
        memcpy(pc+target[i]*width, ps+bestIdx*width, sizeof(T)*width);
    }

    delete []curDist;
    delete []bestDist;
    curDist = bestDist = NULL;
    // printf("done\n");
}

template<class T>
void initClusterCenters(Image<T> &centers, const Image<T> &samples, int clusters,
                        double (*distance)(T*, T*, int))
{
    // printf("using kmeans++ to init cluster centers...\n");
    
    int width = samples.nWidth(), height = samples.nHeight(), idx;
    T *pc, *ps;
    double sum, *dist = new double[height];
    int *target = new int[clusters];
    
    centers.create(width, clusters);
    pc = centers.ptr(), ps = samples.ptr();
    
    // find the first cluster center uniformly randomly
    srand(time(NULL));
    idx = rand() % height;
    memcpy(pc, ps+idx*width, width*sizeof(T));

    // calculate distance
    sum = 0;
    for (int h = 0; h < height; ++h)
    {
        dist[h] = (*distance)(ps+h*width, pc, width);
        sum += dist[h];
    }

    for (int c = 1; c < clusters; ++c)
        target[c-1] = c;
    
    // find remain cluster centers
    findClusterCenters(centers, samples, clusters, dist, sum,
                       target, clusters-1, distance);

    delete []dist;
    delete []target;
    dist = NULL;
    target = NULL;
    
    // printf("done\n");
}

// each row contains one sample
template <class T>
double kmeans(Image<T> &centers, UCImage &labels, const Image<T> &samples,
              int clusters,
              double (*distance)(T*, T*, int)=dist2)
{
    // printf("running kmeans...\n");
    assert(samples.nChannels() == 1 && clusters > 0 && distance != NULL);
    
    int width = samples.nWidth(), height = samples.nHeight(), idx, need2deter;
    const int maxIter = 1000;
    const double threshold = 0.001 * width;
    double best_comp = DBL_MAX, max_shift = DBL_MAX, preshift = DBL_MAX, precomp = DBL_MAX;
    double compactness, shift, dist, tmp, sum;
    double *distArr = new double[height];
    int *counters = new int[clusters];
    int *target = new int[clusters];
    Image<T> nCenters(width, clusters);
    bool empty;
    
    if (centers.isEmpty())
        initClusterCenters(centers, samples, clusters, distance);// kmeans++

    // debug
    // imprint(centers);
    
    labels.create(1, height);
    T *ps = samples.ptr(), *pc = centers.ptr(), *pnc = nCenters.ptr();

    for (int iter = 0; iter < maxIter && max_shift > threshold; ++iter)
    {
        // calc data point to nearest cluster center
        // printf("calc data point to nearest cluster center\n");
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
        empty = false;
        need2deter = 0;
        for (int c = 0; c < clusters; ++c)
        {
            if (counters[c] >= 0) continue;
            empty = true;
            target[need2deter++] = c;
        }

        if (empty)
        {
            printf("find cluster without element, create new cluster...");

            sum = 0;
            for (int h = 0; h < height; ++h)
            {
                distArr[h] = DBL_MAX;
                for (int c = 0; c < clusters; ++c)
                {
                    if (counters[c] <= 0) continue;
                    tmp = (*distance)(ps+h*width, pc+c*width, width);
                    distArr[h] = std::min(distArr[h], tmp);
                }
                sum += distArr[h];
            }
            
            findClusterCenters(centers, samples, clusters, distArr, sum,
                               target, need2deter, distance);
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
            compactness += (*distance)(ps+h*width, pc+labels[h]*width, width);

        if (compactness < best_comp)
        {
            best_comp = compactness;
            nCenters.copyTo(centers);
        }

        // printf("iterations %d, compactness: %.6f, max_shift: %.6f\n", iter, best_comp, max_shift);
        if (fabs(preshift - max_shift) < ESP && fabs(precomp - best_comp) < ESP) break;

        preshift = max_shift;
        precomp = best_comp;
    }

    delete []counters;
    delete []distArr;
    delete []target;
    counters = NULL;
    distArr = NULL;
    target = NULL;
    return best_comp;
}

// adaptive kmeans, try serverl cluster number and get the best one
template <class T>
int kmeans2(Image<T> &centers, UCImage &labels, const Image<T> &samples,
            int start = 2, int end = 20,
            double (*distance)(T*, T*, int)=dist2)
{
    int width = samples.nWidth(), height = samples.nHeight(), clusters, bestClusters;
    double error, preerror = DBL_MAX;
    Image<T> bestCenters;
    UCImage bestLabels(1, height);
    const double logR = log(height);
    const double na = 0.8;
    
    for (clusters = start; clusters <= end; ++clusters)
    {
        centers.release();
        error = kmeans(centers, labels, samples, clusters, distance);
        printf("clusters: %d, compactness: %.6f\n", clusters, error);

        // Schwarz criterion for model selection
        error += na * clusters * width * logR;
        printf("bic: %.6f\n", error);
        if (error < preerror)
        {
            preerror = error;
            bestClusters = clusters;
            centers.copyTo(bestCenters);
            labels.copyTo(bestLabels);
        }
    }

    bestCenters.copyTo(centers);
    bestLabels.copyTo(labels);
    return bestClusters;
}

#endif
