#include "LAOF.h"

void LAOF::EM(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
              DImage &mask1, DImage &mask2,
              const DImage &im1, const DImage &im2, int nIters, int idx)
{
    // init optical flow parameters
    const double as = 0.026;
    const double ap = 0.012;
    const double ratio = 0.75;
    const int minWidth = 20;
    const int nOutIter = 5;//14;
    const int nInIter = 1;
    const int nSORIter = 10;//30;

    int width = im1.nWidth(), height = im1.nHeight();
    DImage tmpu1, tmpv1, tmpu2, tmpv2, flow, rflow, tmp;
    DImage centers, layers1, layers2, fInfo, flowf, rflowf;
    UCImage flowImg;
    int clusters, l, i, iter;//, m1Size, m2Size, mSize;
    // double nas;
    OpticalFlow of;
    MotionLayers ml;

    char buf[256];
    const char *out = "/home/iaml/Projects/exp/lena/out/%d-iter%d-layer%d-%s.jpg";

    layers1.create(width, height);
    layers2.create(width, height);
    mask1.create(width, height);
    mask2.create(width, height);
    flow.create(width, height, 2);
    rflow.create(width, height, 2);
    
    for (clusters = 1, iter = 0; iter < nIters; ++iter)
    {
        printf("EM iterations: %d\n", iter);
        
        // seperate mask
        for (l = 0; l < clusters; ++l)
        {
            // m1Size = 0, m2Size = 0;
            for (i = 0; i < layers1.nElements(); ++i)
            {
                if (fabs(layers1[i]-l) < 0.5)
                {
                    mask1[i] = 1;
                    // m1Size++;
                }
                else mask1[i] = 0;

                if (fabs(layers2[i]-l) < 0.5)
                {
                    mask2[i] = 1;
                    // m2Size++;
                }
                else mask2[i] = 0;
            }
            
            // for test purpose
            sprintf(buf, out, idx, iter, l, "mask1");
            imwrite(buf, mask1);
            sprintf(buf, out, idx, iter, l, "mask2");
            imwrite(buf, mask2);

            // mSize = std::max(m1Size, m2Size);
            // if (clusters == 1) nas = as;
            // else nas = as * log(mSize);
            of.biC2FFlow(tmpu1, tmpv1, tmpu2, tmpv2, im1, im2, mask1, mask2,
                         as, ap, ratio, minWidth, nOutIter, nInIter, nSORIter);

            multiply(tmpu1, mask1);
            multiply(tmpv1, mask1);
            multiply(tmpu2, mask2);
            multiply(tmpv2, mask2);
            
            // sprintf(buf, out, iter, l, "sub-flow");
            // flow2color(flowImg, tmpu1, tmpv1);
            // imwrite(buf, flowImg);
            
            // sprintf(buf, out, iter, l, "sub-rflow");
            // flow2color(flowImg, tmpu2, tmpv2);
            // imwrite(buf, flowImg);

            // merge u, v to flow, ur, vr to rflow
            for (i = 0; i < mask1.nElements(); ++i)
            {
                if (fabs(mask1[i]-1) < ESP)
                {
                    flow[i*2] = tmpu1[i];
                    flow[i*2+1] = tmpv1[i];
                }
                if (fabs(mask2[i]-1) < ESP)
                {
                    rflow[i*2] = tmpu2[i];
                    rflow[i*2+1] = tmpv2[i];
                }
            }
        }

        // for test purpose
        sprintf(buf, out, idx, iter, 0, "flow");
        flow2color(flowImg, flow);
        imwrite(buf, flowImg);
            
        sprintf(buf, out, idx, iter, 0, "flowr");
        flow2color(flowImg, rflow);
        imwrite(buf, flowImg);
        
        // clustering
        ml.flowInfo(flowf, flow);
        ml.flowInfo(rflowf, rflow);

        if (centers.isEmpty())
        {
            printf("first time to run cluster, automatically determine # clusters\n");
            clusters = ml.cluster(centers, layers1, flowf, width, height, 2, 5, true);
        } else {
            createCenterByLabel(centers, clusters, layers1, flowf);
            ml.cluster(centers, layers1, flowf, width, height, clusters, clusters);
        }

        sprintf(buf, out, idx, iter, 0, "layers1");
        flow2color(flowImg, layers1, layers1);
        imwrite(buf, flowImg);

        ml.coverLabels(tmp, im1, flowImg);
        sprintf(buf, out, idx, iter, 0, "merge-layers1");
        imwrite(buf, tmp);

        // refine layers
        ml.refine(layers1, clusters, im1, flow, centers, flowf);

        sprintf(buf, out, idx, iter, 0, "refine-layers1");
        flow2color(flowImg, layers1, layers1);
        imwrite(buf, flowImg);

        ml.coverLabels(tmp, im1, flowImg);
        sprintf(buf, out, idx, iter, 0, "merge-layers1-refine");
        imwrite(buf, tmp);

        // transfer labels from im1 to im2
        transferLabel(layers2, layers1, flow);
        createCenterByLabel(centers, clusters, layers2, rflowf);
        ml.cluster(centers, layers2, rflowf, width, height, clusters, clusters);

        sprintf(buf, out, idx, iter, 0, "layers2");
        flow2color(flowImg, layers2, layers2);
        imwrite(buf, flowImg);

        ml.coverLabels(tmp, im2, flowImg);
        sprintf(buf, out, idx, iter, 0, "merge-layers2");
        imwrite(buf, tmp);

        // refine layers
        ml.refine(layers2, clusters, im2, rflow, centers, rflowf);

        sprintf(buf, out, idx, iter, 0, "refine-layers2");
        flow2color(flowImg, layers2, layers2);
        imwrite(buf, flowImg);

        ml.coverLabels(tmp, im2, flowImg);
        sprintf(buf, out, idx, iter, 0, "merge-layers2-refine");
        imwrite(buf, tmp);
    }

    u1.create(width, height), v1.create(width, height);
    u2.create(width, height), v2.create(width, height);
    for (int i = 0; i < flow.nSize(); ++i)
    {
        u1[i] = flow[i*2];
        v1[i] = flow[i*2+1];
        u2[i] = rflow[i*2];
        v2[i] = rflow[i*2+1];
    }
   
    printf("EM done!\n");
}

void LAOF::createCenterByLabel(DImage &centers, int clusters,
                               const DImage &labels, const DImage &samples)
{
    // printf("create center by label... ");
    
    int *count = new int[clusters];
    int width = samples.nWidth(), lsize = labels.nSize(), idx;
    
    centers.create(width, clusters);
    memset(count, 0, sizeof(int)*clusters);

    for (int i = 0; i < lsize; ++i)
    {
        idx = labels[i];
        count[idx]++;
        for (int w = 0; w < width; ++w)
            centers[idx*width+w] += samples[i*width+w];
    }

    for (int i = 0; i < clusters; ++i)
        for (int w = 0; w < width; ++w)
            centers[i*width+w] /= (double)count[i];
    
    delete []count;
    // printf("done\n");
}

void LAOF::transferLabel(DImage &dst, const DImage &src, const DImage &flow)
{
    // printf("transfer label... ");

    int width = src.nWidth(), height = src.nHeight(), channels = src.nChannels();
    int nx, ny, offset, noffset;

    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            nx = w + (flow[offset*2] + 0.5);
            ny = h + (flow[offset*2+1] + 0.5);
            if (nx > width || nx < 0 || ny > height || ny < 0)
                continue;

            offset *= channels;
            noffset = (ny * width + nx) * channels;
            for (int k = 0; k < channels; ++k)
                dst[noffset+k] = src[offset+k];
        }
    }

    // printf("done\n");
}
