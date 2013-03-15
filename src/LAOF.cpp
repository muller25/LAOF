#include "LAOF.h"
#include "OpticalFlow.h"
#include "MotionLayers.h"
#include "ImageIO.h"
#include "Flow2Color.h"
#include "OFPara.h"

// init optical flow parameters
const double as = 0.026;
const double ap = 0.012;
const double ratio = 0.75;
const int minWidth = 20;
const int nBiIter = 5;//14;
const int nIRLSIter = 1;
const int nSORIter = 10;//30;

const char *out = "/home/iaml/Projects/exp/lena/out/%d-iter%d-layer%d-%s.jpg";
int frameID = 0;

void LAOF::EM(std::vector<DImage> &u, std::vector<DImage> &v,
              std::vector<DImage> &ur, std::vector<DImage> &vr,
              std::vector<DImage> &masks, const std::vector<DImage> &im,
              int curIdx, int nIters)
{
    assert(u.size() == v.size() && v.size() == ur.size() && ur.size() == masks.size() &&
           masks.size() == im.size() && im.size() == 3);

    int prev = (curIdx+2) % 3, next = (curIdx+1) % 3;
    DImage prfInfo, emptyfInfo, fInfo, rfInfo, flow, rflow, tmp;
    DImage centers1, centers2;
    UCImage ucimg;
    std::vector<DImage> vec;
    std::vector<OFPara> paras;
    MotionLayers ml;
    OpticalFlow of;
    int width = im[0].nWidth(), height = im[0].nHeight(), numOfLabels = 1;
    char buf[256];
    
    // init mask
    printf("******** laof: %d->%d->%d ********\n", prev, curIdx, next);
    masks[curIdx].create(width, height);
    masks[next].create(width, height);
    
    // check if previous flow info is available
    emptyfInfo.create(2, width*height);
    if (!ur[prev].isEmpty() && !vr[prev].isEmpty())
    {
        printf("previous flow info is available\n");
        vec.push_back(ur[prev]);
        vec.push_back(vr[prev]);
        mergec(tmp, vec);
        ml.flowInfo(prfInfo, tmp);
    } else {
        printf("no previous flow info is available\n");
        prfInfo.create(2, width*height);
    }

    // start EM iterations
    for (int iter = 0; iter < nIters; ++iter)
    {
        printf("EM iterations: %d\n", iter);

        // split images, masks into pieces
        OFPara::split(paras, numOfLabels, im[curIdx], im[next],
                      masks[curIdx], masks[next]);

        // loop over layers and run biDir optical flow
        printf("loop over %d layers and run biDir optical flow\n", numOfLabels);
        for (int l = 0; l < numOfLabels; ++l)
        {
            // for test only
            sprintf(buf, out, frameID, iter, l, "mask1");
            imwrite(buf, paras[l].mask1);
            
            sprintf(buf, out, frameID, iter, l, "mask2");
            imwrite(buf, paras[l].mask2);

            // sprintf(buf, out, frameID, iter, l, "im1");
            // imwrite(buf, paras[l].im1);

            // sprintf(buf, out, frameID, iter, l, "im2");
            // imwrite(buf, paras[l].im2);

            of.biC2FFlow(paras[l].u1, paras[l].v1, paras[l].u2, paras[l].v2,
                         paras[l].im1, paras[l].im2,
                         paras[l].mask1, paras[l].mask2,
                         as, ap, ratio, minWidth, nBiIter, nIRLSIter, nSORIter);

            // for test only
            printf("****** sub flow ******\n");
            sprintf(buf, out, frameID, iter, l, "sub-flow");
            flow2color(ucimg, paras[l].u1, paras[l].v1);
            imwrite(buf, ucimg);
            sprintf(buf, out, frameID, iter, l, "sub-rflow");
            flow2color(ucimg, paras[l].u2, paras[l].v2);
            imwrite(buf, ucimg);

            warpImage(tmp, paras[l].im1, paras[l].im2, paras[l].u1, paras[l].v1);
            sprintf(buf, out, frameID, iter, l, "sub-warp");
            imwrite(buf, tmp);

            warpImage(tmp, paras[l].im2, paras[l].im1, paras[l].u2, paras[l].v2);
            sprintf(buf, out, frameID, iter, l, "sub-rwarp");
            imwrite(buf, tmp);
        }

        OFPara::restore(flow, rflow, paras, width, height);

        // for test only
        printf("****** global flow ******\n");
        sprintf(buf, out, frameID, iter, 0, "flow");
        flow2color(ucimg, flow);
        imwrite(buf, ucimg);
        sprintf(buf, out, frameID, iter, 0, "rflow");
        flow2color(ucimg, rflow);
        imwrite(buf, ucimg);

        // extract motion layers
        printf("****** extract motion layer ******\n");

        ml.flowInfo(fInfo, flow);
        vec.clear();
        vec.push_back(fInfo);
        vec.push_back(prfInfo);
        mergew(fInfo, vec);
        
        // kmeans
        if (numOfLabels == 1)
            numOfLabels = ml.cluster(centers1, masks[curIdx], fInfo,
                                     width, height, 4, 4, 15, true);
        else
            ml.cluster(centers1, masks[curIdx], fInfo, width, height,
                       numOfLabels, numOfLabels);

        // for test only
        sprintf(buf, out, frameID, iter, 0, "layers1");
        flow2color(ucimg, masks[curIdx], masks[curIdx]);
        imwrite(buf, ucimg);

        ml.coverLabels(tmp, im[curIdx], ucimg);
        sprintf(buf, out, frameID, iter, 0, "merge-layers1");
        imwrite(buf, tmp);

        // refine layers of im1
        ml.refine(masks[curIdx], numOfLabels, im[curIdx], flow, centers1, fInfo);
        ml.createCenterByLabels(centers1, numOfLabels, masks[curIdx], fInfo);
        
        // for test only
        sprintf(buf, out, frameID, iter, 0, "layers1-refine");
        flow2color(ucimg, masks[curIdx], masks[curIdx]);
        imwrite(buf, ucimg);

        ml.coverLabels(tmp, im[curIdx], ucimg);
        sprintf(buf, out, frameID, iter, 0, "merge-layers1-refine");
        imwrite(buf, tmp);

        // transfer labels from im1 to im2
        ml.flowInfo(rfInfo, rflow);
        vec.clear();
        vec.push_back(rfInfo);
        vec.push_back(emptyfInfo);
        mergew(rfInfo, vec);
        
        transferLabels(masks[next], masks[curIdx], flow);
        ml.createCenterByLabels(centers2, numOfLabels, masks[next], rfInfo);
        ml.cluster(centers2, masks[next], rfInfo, width, height,
                   numOfLabels, numOfLabels);

        // for test only
        sprintf(buf, out, frameID, iter, 0, "layers2");
        flow2color(ucimg, masks[next], masks[next]);
        imwrite(buf, ucimg);

        ml.coverLabels(tmp, im[next], ucimg);
        sprintf(buf, out, frameID, iter, 0, "merge-layers2");
        imwrite(buf, tmp);

        // refine layers of im2
        ml.refine(masks[next], numOfLabels, im[next], rflow, centers2, rfInfo);
        ml.createCenterByLabels(centers2, numOfLabels, masks[next], rfInfo);
        
        // for test only
        sprintf(buf, out, frameID, iter, 0, "layers2-refine");
        flow2color(ucimg, masks[next], masks[next]);
        imwrite(buf, ucimg);

        ml.coverLabels(tmp, im[next], ucimg);
        sprintf(buf, out, frameID, iter, 0, "merge-layers2-refine");
        imwrite(buf, tmp);
    }
    // end of EM iterations

    // set return parameters
    int offset;
    u[curIdx].create(width, height);
    v[curIdx].create(width, height);
    ur[curIdx].create(width, height);
    vr[curIdx].create(width, height);
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = h * width + w;
            u[curIdx][offset] = flow[offset*2];
            v[curIdx][offset] = flow[offset*2+1];
            ur[curIdx][offset] = rflow[offset*2];
            vr[curIdx][offset] = rflow[offset*2+1];
        }
    }

    frameID++;
}

void LAOF::transferLabels(DImage &dst, const DImage &srcLabels, const DImage &flow)
{
    // printf("transfer label... ");

    int width = srcLabels.nWidth(), height = srcLabels.nHeight(), channels = srcLabels.nChannels();
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
                dst[noffset+k] = srcLabels[offset+k];
        }
    }

    // printf("done\n");
}
