#ifndef _ML_H
#define _ML_H

#include <vector>
#include <ctime>
using std::vector;

#define ESP 1e-6

template <class T>
inline double dist2(const T *p1, const T *p2, int start, int end)
{
    double d = 0;
    for (int i = start; i < end; ++i)
        d += (p1[i] - p2[i]) * (p1[i] - p2[i]);

    return sqrt(d);
}

template <class T>
int kmeans2(vector< vector<T> > &centers, vector<int> &labels,
            const vector< vector<T> > &samples,
            int start = 2, int end = 20, double na = 2,
            double (*distance)(const T*, const T*, int, int)=dist2)
{
    int ndata = samples[0].size(), nsample = samples.size();
    int bestClusters = -1;
    double error, preerror = DBL_MAX;
    vector<T> bestCenters;
    vector<int> bestLabels;
    const double logR = log(nsample);
    
    for (int ncluster = start; ncluster <= end; ++ncluster)
    {
        error = kmeans(centers, labels, samples, ncluster, distance);
        // printf("clusters: %d, compactness: %.6f\n", ncluster, error);

        // Schwarz criterion for model selection
        error += na * ncluster * ndata * logR;
        // printf("bic: %.6f\n", error);
        if (error < preerror)
        {
            preerror = error;
            bestClusters = ncluster;
            bestCenters = centers;
            bestLabels = labels;
        }
    }

    centers = bestCenters;
    labels = bestLabels;
    return bestClusters;
}

template <class T>
double kmeans(vector< vector<T> > &centers, vector<int> &labels,
              const vector< vector<T> > &samples,
              int ncluster,
              double (*distance)(const T*, const T*, int, int)=dist2)
{
    assert(ncluster > 0 && distance != NULL);

    int nsample = samples.size(), ndata = samples[0].size(), idx, need2deter;
    const int maxIter = 1000;
    const double threshold = 0.001 * ndata;
    double best_comp = DBL_MAX, precomp = DBL_MAX;
    double max_shift = DBL_MAX, preshift = DBL_MAX;
    double compactness, shift, dist, tmp, sum;
    vector<double> distArr(nsample, 0);
    vector<int> counters(ncluster), target(ncluster, 0);
    vector< vector<T> > nCenters;
    bool empty;
    
    if ((int)centers.size() != ncluster)
        initClusterCenters(centers, samples, ncluster, distance);// kmeans++

    labels.resize(nsample);
    for (int iter = 0; iter < maxIter && max_shift > threshold; ++iter)
    {
        // calc data point to nearest cluster center
        nCenters.resize(ncluster);
        for (int i = 0; i < ncluster; ++i)
        {
            nCenters[i].resize(ndata);
            memset(nCenters[i].data(), 0, sizeof(T) * ndata);
        }
        
        memset(counters.data(), 0, sizeof(int) * ncluster);
        for (int i = 0; i < nsample; ++i)
        {
            dist = DBL_MAX, idx = -1;
            for (int c = 0; c < ncluster; ++c)
            {
                tmp = (*distance)(samples[i].data(), centers[c].data(), 0, ndata);
                if (tmp < dist)
                {
                    dist = tmp;
                    idx = c;
                }
            }

            // assign label
            labels[i] = idx;
            counters[idx]++;
            for (int w = 0; w < ndata; ++w)
                nCenters[idx][w] += samples[i][w];
        }

        // check whether there is cluster without element
        empty = false;
        need2deter = 0;
        for (int c = 0; c < ncluster; ++c)
        {
            if (counters[c] > 0) continue;
            empty = true;
            target[need2deter++] = c;
        }

        if (empty)
        {
            // printf("find cluster without element, create new cluster...");

            sum = 0;
            for (int h = 0; h < nsample; ++h)
            {
                distArr[h] = DBL_MAX;
                for (int c = 0; c < ncluster; ++c)
                {
                    if (counters[c] <= 0) continue;
                    tmp = (*distance)(samples[h].data(), centers[c].data(), 0, ndata);
                    distArr[h] = std::min(distArr[h], tmp);
                }
                sum += distArr[h];
            }
            
            findClusterCenters(centers, samples, ncluster, distArr, sum,
                               target, need2deter, distance);
            // printf("done\n");
        }
        
        // calc new cluster centers && center shift
        max_shift = 0;
        for (int c = 0; c < ncluster; ++c)
        {
            shift = 0;
            for (int w = 0; w < ndata; ++w)
            {
                nCenters[c][w] /= counters[c];
                shift += fabs(nCenters[c][w] - centers[c][w]);
            }

            if (max_shift < shift) max_shift = shift;
        }

        // calc compactness of clusters
        compactness = 0;
        for (int h = 0; h < nsample; ++h)
            compactness += (*distance)(samples[h].data(), centers[labels[h]].data(), 0, ndata);

        if (compactness < best_comp)
        {
            best_comp = compactness;
            nCenters = centers;
        }

        // printf("iterations %d, compactness: %.6f, max_shift: %.6f\n", iter, best_comp, max_shift);
        if (fabs(preshift - max_shift) < ESP && fabs(precomp - best_comp) < ESP) break;

        preshift = max_shift;
        precomp = best_comp;
    }

    return best_comp;
}

template<class T>
void findClusterCenters(vector< vector<T> > &centers,
                        const vector< vector<T> > &samples, int ncluster,
                        vector<double> &dist, double sum,
                        const vector<int> &target, int need2deter,
                        double (*distance)(const T*, const T*, int, int))
{
    const int trials = 5;
    int ndata = samples[0].size(), nsample = samples.size(), curIdx, bestIdx, h;
    double curSum, bestSum, tmp;
    vector<double> curDist(nsample), bestDist(nsample);
        
    srand(time(NULL));
    for (int i = 0; i < need2deter; ++i)
    {
        // try to find cluster centers
        bestSum = DBL_MAX, bestIdx = -1;
        for (int t = 0; t < trials; ++t)
        {
            // proportional random
            curSum = sum * (double)rand() / RAND_MAX;
            // printf("sum: %.6f, curSum %.6f\n", sum, curSum);
            
            for (curIdx = 0; curIdx < nsample; ++curIdx)
                if ((curSum -= dist[curIdx]) <= 0) break;

            // printf("curSum %.6f, curIdx %d\n", curSum, curIdx);
            assert(curIdx < nsample);

            // find minimum distance to cluster centers
            for (curSum = 0, h = 0; h < nsample; ++h)
            {
                tmp = (*distance)(samples[h].data(), samples[curIdx].data(), 0, ndata);
                curDist[h] = std::min(dist[h], tmp);
                curSum += curDist[h];
            }

            // is this one better than previous one?
            if (curSum < bestSum)
            {
                bestSum = curSum, bestIdx = curIdx;
                bestDist = curDist;
            }
        }
        
        assert(bestIdx >= 0);
        sum = bestSum;
        dist = bestDist;
        centers[target[i]] = samples[bestIdx];
    }
}

template<class T>
void initClusterCenters(vector< vector<T> > &centers,
                        const vector< vector<T> > &samples, int ncluster,
                        double (*distance)(const T*, const T*, int, int))
{
    int ndata = samples[0].size(), nsample = samples.size(), idx;

    centers.resize(ncluster);
    for (int c = 0; c < ncluster; ++c)
    {
        centers[c].resize(ndata);
        memset(centers[c].data(), 0, sizeof(T)*ndata);
    }
    
    // find the first cluster center uniformly randomly
    srand(time(NULL));
    idx = rand() % nsample;
    centers[0] = samples[idx];

    // calculate distance
    vector<double> dist(nsample);
    double sum = 0;
    for (int h = 0; h < nsample; ++h)
    {
        dist[h] = (*distance)(samples[h].data(), centers[0].data(), 0, ndata);
        sum += dist[h];
    }

    assert(!isnan(sum));
    
    vector<int> target(nsample);
    for (int c = 1; c < ncluster; ++c)
        target[c-1] = c;
    
    // find remain cluster centers
    findClusterCenters(centers, samples, ncluster, dist, sum,
                       target, ncluster-1, distance);
}

#endif
