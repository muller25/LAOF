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
            data[i*labels+l] = 2 * mydist(features.ptr()+i*cwidth, centers.ptr()+l*cwidth, 0, cwidth);
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
    const double penalty = 3;
    const double sigma = 3;
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

void MotionLayers::scluster(DImage &centers, DImage &layers, int nlabels,
                            const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const int wsize = 1;
    const double threshold = ESP;
    const int neighbour = 5 * wsize;

    int width = u.nWidth(), height = u.nHeight(), offset, h, w;
    
    DImage covUV;
    int trust = 0, totalPoints = u.nSize();
    int *t2id = new int[totalPoints];
    int *id2t = new int[totalPoints];
    
    for (int i = 0; i < totalPoints; ++i)
        id2t[i] = -1;
    
    // count trusted points    
    covariance(covUV, u, v, wsize);
    for (h = wsize; h < height; h += wsize)
    {
        for (w = wsize; w < width; w += wsize)
        {
            offset = h * width + w;
            if (fabs(covUV[offset]) <= threshold)
            {
                t2id[trust] = offset;
                id2t[offset] = trust++;
            }
        }
    }

    printf("trusted points: %d\n", trust);

    // construct similarity matrix
    printf("construct similarity matrix... ");

    int nh, nw, noffset, j;
    double cost, tmpu, tmpv;
    DImage graph(trust, trust);
    for (int i = 0; i < trust; ++i)
    {
        h = t2id[i] / width;
        w = t2id[i] % width;

        // search for trusted neighbours
        for (int hh = -neighbour; hh <= neighbour; hh += wsize)
        {
            nh = h + hh;
            if (nh < 0 || nh >= height) continue;

            for (int ww = -neighbour; ww <= neighbour; ww += wsize)
            {
                nw = w + ww;
                if (nw < 0 || nw >= width) continue;

                noffset = nh * width + nw;
                if (id2t[noffset] == -1) continue;

                j = id2t[noffset];
                offset = t2id[i];
                tmpu = u[offset] - u[noffset];
                tmpv = v[offset] - v[noffset];
                cost = sqrt(tmpu*tmpu + tmpv*tmpv);
                graph[i*trust+j] = cost;
            }
        }
    }

    printf("done\n");
    delete []id2t;

    // spectral cluster
    printf("running spectral cluster...");

    UCImage ucimg;
    SpectralCluster(ucimg, graph, nlabels);

    printf("done\n");
    
    int *count = new int[nlabels];
    int label;
    
    memset(count, 0, sizeof(int)*nlabels);
    centers.create(2, nlabels);
    layers.create(width, height, 1, -1);
    for (int i = 0; i < trust; ++i)
    {
        label = ucimg[i];
        offset = t2id[i];
        layers[offset] = label;
        count[label]++;
        centers[label*2] += u[offset];
        centers[label*2+1] += v[offset];
    }

    for (int i = 0; i < nlabels; ++i)
    {
        centers[i*2] /= count[i];
        centers[i*2+1] /= count[i];
    }
    
    delete []count;
    delete []t2id;
}

void MotionLayers::refine(DImage &centers, DImage &layers, int nlabels,
                          const DImage &im1, const DImage &im2,
                          const DImage &u, const DImage &v)
{
    assert(centers.ptr() != NULL);
    
    int size = im1.nSize();
    int width = im1.nWidth(), height = im1.nHeight(), channels = im1.nChannels();
    DImage warp, features(2, size), extra(width, height, channels*2);

    warpImage(warp, im1, im2, u, v);
    for (int i = 0; i < size; ++i)
    {
        features[i*2] = u[i];
        features[i*2+1] = v[i];

        for (int k = 0; k < channels; ++k)
            extra[i*channels*2 + k] = im1[i*channels + k];

        for (int k = 0; k < channels; ++k)
            extra[i*channels*2+channels+k] = warp[i*channels+k];
    }

    double *data = NULL;    
    GCoptimization *gc = NULL;
    try{
		gc = new GCoptimizationGridGraph(width, height, nlabels);
        data = new double[nlabels*size];
        gc->setSmoothCost(MotionLayers::smoothFn, extra.ptr());

        for (int iter = 0; iter < 1; ++iter)
        {
            dataFn(data, nlabels, centers, features);
            gc->setDataCost(data);
        
            printf("Before optimization energy is %.6f\n",gc->compute_energy());
            gc->expansion(2);
            // gc->swap(2);
            printf("After optimization energy is %.6f\n",gc->compute_energy());

            for (int i = 0; i < size; ++i)
                layers[i] = gc->whatLabel(i);

            createCenterByLabels(centers, nlabels, layers, features);
        }
	}
	catch (GCException e){
		e.Report();
	}

    if (gc != NULL) delete gc;
    if (data != NULL) delete []data;
}

void MotionLayers::dataFn(double *data, int nlabels,
                          const DImage &centers, const DImage &features)
{
    int size = features.nHeight(), cwidth = centers.nWidth();
    double *pc = centers.ptr(), *pf = features.ptr();
    
    for (int i = 0; i < size; ++i)
        for (int l = 0; l < nlabels; ++l)
            data[i*nlabels+l] = dist2(pf+i*cwidth, pc+l*cwidth, 0, cwidth);
}
