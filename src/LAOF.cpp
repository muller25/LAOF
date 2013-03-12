#include "LAOF.h"

void LAOF::EM(DImage &u1, DImage &v1, DImage &u2, DImage &v2,
              DImage &mask1, DImage &mask2,
              const DImage &im1, const DImage &im2, int nIters)
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
    DImage tmpu1, tmpv1, tmpu2, tmpv2;
    DImage centers, layers1, layers2, flow, rflow, fInfo, flowf, rflowf;
    UCImage flowImg;
    int clusters, l, i, iter;
    OpticalFlow of;
    MotionLayers ml;

    char buf[256];
    const char *out = "/home/iaml/Projects/exp/lena/out/iter%d-layer%d-%s.jpg";
    const char *outWarp = "/home/iaml/Projects/exp/lena/out/%swarp%03d.jpg";
    
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
            for (i = 0; i < layers1.nElements(); ++i)
            {
                if (fabs(layers1[i]-l) < 0.5) mask1[i] = 1;
                else mask1[i] = 0;

                if (fabs(layers2[i]-l) < 0.5) mask2[i] = 1;
                else mask2[i] = 0;
            }
            
            // for test purpose
            sprintf(buf, out, iter, l, "mask1");
            imwrite(buf, mask1);
            sprintf(buf, out, iter, l, "mask2");
            imwrite(buf, mask2);

            of.biC2FFlow(tmpu1, tmpv1, tmpu2, tmpv2, im1, im2, mask1, mask2,
                         as, ap, ratio, minWidth, nOutIter, nInIter, nSORIter);

            multiply(tmpu1, mask1);
            multiply(tmpv1, mask1);
            multiply(tmpu2, mask2);
            multiply(tmpv2, mask2);
            
            sprintf(buf, out, iter, l, "sub-flow");
            flow2color(flowImg, tmpu1, tmpv1);
            imwrite(buf, flowImg);
            
            sprintf(buf, out, iter, l, "sub-rflow");
            flow2color(flowImg, tmpu2, tmpv2);
            imwrite(buf, flowImg);

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
        sprintf(buf, out, iter, iter, "flow");
        flow2color(flowImg, flow);
        imwrite(buf, flowImg);
            
        sprintf(buf, out, iter, iter, "rflow");
        flow2color(flowImg, rflow);
        imwrite(buf, flowImg);
        
        // clustering
        ml.flowInfo(flowf, flow);
        ml.flowInfo(rflowf, rflow);

        for (int j = 0; j < 2; ++j)
        {
            printf("clustering iteration: %d\n", j);

            if (centers.isEmpty())
            {
                printf("first time to run cluster, automatically determine # clusters\n");
                clusters = ml.cluster(centers, layers1, flowf, width, height, 2, 5);
            }
            else
                ml.cluster(centers, layers1, flowf, width, height, clusters, clusters);
            
            transferLabel(layers2, layers1, flow);
            createCenterByLabel(centers, clusters, layers2, rflowf);

            sprintf(buf, out, iter, j, "layers1");
            flow2color(flowImg, layers1, layers1);
            imwrite(buf, flowImg);

            sprintf(buf, out, iter, j, "layers2t");
            flow2color(flowImg, layers2, layers2);
            imwrite(buf, flowImg);
        
            ml.cluster(centers, layers2, rflowf, width, height, clusters, clusters);
            transferLabel(layers1, layers2, rflow);
            createCenterByLabel(centers, clusters, layers1, flowf);

            sprintf(buf, out, iter, j, "layers2");
            flow2color(flowImg, layers2, layers2);
            imwrite(buf, flowImg);

            sprintf(buf, out, iter, j, "layers1t");
            flow2color(flowImg, layers1, layers1);
            imwrite(buf, flowImg);
        }
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

    DImage warp;
    of.warpImage(warp, im1, im2, u1, v1);
    sprintf(buf, outWarp, "", 0);
    imwrite(buf, warp);

    of.warpImage(warp, im2, im1, u2, v2);
    sprintf(buf, outWarp, "r", 1);
    imwrite(buf, warp);
   
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
