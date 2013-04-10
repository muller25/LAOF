#include "GCoptimization.h"
#include "MotionLayers.h"
#include "ML.h"
#include "Maths.h"
#include "Utils.h"

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
                       MotionLayers::kmdist);

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
                       MotionLayers::kmdist);

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
            data[i*labels+l] = 2 * kmdist(features.ptr()+i*cwidth, centers.ptr()+l*cwidth, 0, cwidth);
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

// kmeans cluster, -1 for untrusted areas
void MotionLayers::kcluster(DImage &centers, DImage &layers, int nlabels,
                            const DImage &im1, const DImage &im2,
                            const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const int wsize = 2;
    const double threshold = ESP;
    int height = u.nHeight(), width = u.nWidth();
    
    DImage covUV;
    int trust = 0, totalPoints = u.nSize();
    int *t2id = new int[totalPoints];
    
    // count trusted points    
    covariance(covUV, u, v, wsize);
    for (int i = 0; i < totalPoints; ++i)
        if (covUV[i] <= threshold)
            t2id[trust++] = i;
    
    printf("trusted points: %d\n", trust);

    // generate trusted samples
    int swidth = 4, offset;
    DImage samples(swidth, trust);
    for (int i = 0; i < trust; ++i)
    {
        offset = t2id[i];
        samples[i*swidth] = offset % width;
        samples[i*swidth+1] = offset / width;
        samples[i*swidth+2] = u[offset];
        samples[i*swidth+3] = v[offset];
    }

    // run kmeans, data contains spatial info and motion info
    UCImage ucimg;
    kmeans(centers, ucimg, samples, nlabels, MotionLayers::kmdist);

    // -1 for untrusted point
    layers.create(width, height, 1, -1);
    for (int i = 0; i < trust; ++i)
        layers[t2id[i]] = ucimg[i];
    
    delete []t2id;
}

// spectral cluster, very slow
void MotionLayers::scluster(DImage &centers, DImage &layers, int nlabels,
                            const DImage &u, const DImage &v)
{
    assert(u.match3D(v) && u.nChannels() == 1);

    const int wsize = 2;
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

// distance for kmeans
double MotionLayers::kmdist(double *p1, double *p2, int start, int end)
{
    double cost = 0, dist = 0;
    int idx = 0;
    
    // spatial info
    dist = dist2(p1, p2, idx, idx+sWidth) + 9;
    idx += sWidth;
    
    // flow info
    cost += log10(dist) * dist2(p1, p2, idx, idx+fWidth);
    idx += fWidth;
    
    return cost;
}

// graph cut refine
void MotionLayers::refine(DImage &centers, DImage &layers, int nlabels,
                          const DImage &im1, const DImage &im2,
                          const DImage &u, const DImage &v)
{
    assert(centers.ptr() != NULL);

    int size = im1.nSize(), width = im1.nWidth(), height = im1.nHeight();
    DImage im1lab, im2lab;
    double *data = NULL;    
    GCoptimization *gc = NULL;

    // color space convertion
    BGR2Lab(im1lab, im1);
    BGR2Lab(im2lab, im2);
    
    try{
		gc = new GCoptimizationGridGraph(width, height, nlabels);

        // init data term
        data = new double[nlabels*size];
        dataFn(data, nlabels, centers, layers, u, v, im1lab, im2lab);
        gc->setDataCost(data);
        gc->setSmoothCost(MotionLayers::smoothFn, im1lab.ptr());
        
        printf("Before optimization energy is %.6f\n",gc->compute_energy());
        gc->expansion(3);
        // gc->swap(3);
        printf("After  optimization energy is %.6f\n",gc->compute_energy());

        printf("smooth: %.6f .. %.6f\n", mins, maxs);
        
        for (int i = 0; i < size; ++i)
            layers[i] = gc->whatLabel(i);
	}
	catch (GCException e){
		e.Report();
	}

    if (gc != NULL) delete gc;
    if (data != NULL) delete []data;
}

// D(lp) = -w1 * ln(Cp | lp) - w2 * ln(fp | lp)
void MotionLayers::dataFn(double *data, int nlabels, const DImage &centers, const DImage &layers,
                          const DImage &u, const DImage &v, const DImage &im1, const DImage &im2)
{
    printf("calculating data term...\n");

    const double omegap = 1;
    const double omegac = 0.5;
    const double omegat = 0;
    const double omegaf = 0.085;
    double mind = DBL_MAX;
    double maxd = 0;
    
    int size = im1.nSize(), width = im1.nWidth(), channels = im1.nChannels(), offset;
    DImage labProb, omProb, om, warp;
    UCImage mask;

    warpImage(warp, im1, im2, u, v);
    substract(warp, im1);
    
    // 根据u，v计算方向-强度(om)图，并归一化强度
    double maxrad = -1;
    om.create(u.nWidth(), u.nHeight(), 2);
    for (int i = 0; i < size; ++i)
    {
        offset = i * 2;

        // 方向, [-PI, PI]
        om[offset] = atan2(v[i], u[i]);

        // 强度， [0, 1]
        om[offset+1] = sqrt(u[i] * u[i] + v[i] * v[i]);
        maxrad = std::max(om[offset+1], maxrad);
    }
    if (fabs(maxrad) > ESP)
        for (int i = 0; i < size; ++i)
            om[i*2+1] /= maxrad;

    // 计算data term
    double lp, omp, cx, cy, x, y, pCost, tCost;
    int cwidth = centers.nWidth();
    for (int l = 0; l < nlabels; ++l)
    {
        genLayerMask(mask, layers, l);
        LabComfirmity(labProb, im1, mask);
        OMComfirmity(omProb, om, mask);
        cx = centers[l*cwidth];
        cy = centers[l*cwidth+1];
        for (int i = 0; i < size; ++i)
        {
            // truncate probability smaller ESP 
            if (fabs(labProb[i]) < ESP) lp = -log(ESP);
            else lp = -log(labProb[i]);

            // truncate probability smaller ESP
            if (fabs(omProb[i]) < ESP) omp = -log(ESP);
            else omp = -log(omProb[i]);

            // position
            x = i % width - cx;
            y = i / width - cy;
            pCost = log10(sqrt(x * x + y * y) + 10);

            // motion
            tCost = 0;
            for (int c = 0; c < channels; ++c)
                tCost = warp[i*channels+c] * warp[i*channels+c];
            tCost = sqrt(tCost);
            
            data[i*nlabels+l] = omegap * pCost + omegac * lp + omegaf * omp + omegat * tCost;
            
            mind = std::min(mind, data[i*nlabels+l]);
            maxd = std::max(maxd, data[i*nlabels+l]);
        }
    }

    printf("data: %.6f .. %.6f\n", mind, maxd);
    printf("done\n");
}

double MotionLayers::mins = DBL_MAX;
double MotionLayers::maxs = 0;

double MotionLayers::smoothFn(int p1, int p2, int l1, int l2, void *pData)
{
    const double weight = 8;
    const double tao = 10;
    double *ptr = (double *)pData;
    double cost;
    
    if (l1 == l2) return 0;

    cost = dist2(ptr+p1*3, ptr+p2*3, 0, 3);
    cost = PI / 2 - atan(cost - tao);
    cost = weight * cost;

    mins = std::min(mins, cost);
    maxs = std::max(maxs, cost);
    
    return cost;
}
