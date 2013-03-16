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

int frameID = 0;
char buf[256];
char LAOF::out[256];

void LAOF::EM(std::vector<DImage> &u, std::vector<DImage> &v,
              std::vector<DImage> &ur, std::vector<DImage> &vr,
              std::vector<DImage> &masks, const std::vector<DImage> &im,
              int curIdx, int nIters)
{
    assert(u.size() == v.size() && v.size() == ur.size() && ur.size() == masks.size() &&
           masks.size() == im.size() && im.size() == 3);

    // int prev = (curIdx+2) % 3;

    int next = (curIdx+1) % 3;
    int width = im[0].nWidth(), height = im[0].nHeight(), numOfLabels = 1;
    DImage tmp, centers1, centers2, flow, rflow;
    UCImage ucimg;
    std::vector<OFPara> paras;
    OpticalFlow of;
    
    // init mask
    masks[curIdx].create(width, height);
    masks[next].create(width, height);
    centers2.create(4, width*height);
    
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
            if (l > 0)
            {
                sprintf(buf, out, frameID, iter, l, "mask1");
                imwrite(buf, paras[l].mask1);
            
                sprintf(buf, out, frameID, iter, l, "mask2");
                imwrite(buf, paras[l].mask2);
            }
            
            of.biC2FFlow(paras[l].u1, paras[l].v1, paras[l].u2, paras[l].v2,
                         paras[l].im1, paras[l].im2,
                         paras[l].mask1, paras[l].mask2,
                         as, ap, ratio, minWidth, nBiIter, nIRLSIter, nSORIter);

            // for test only
            printf("****** sub flow ******\n");
            printf("u1 %.6f .. %.6f, v1 %.6f .. %.6f\n", paras[l].u1.min(), paras[l].u1.max(), paras[l].v1.min(), paras[l].v1.max());
            printf("u2 %.6f .. %.6f, v2 %.6f .. %.6f\n", paras[l].u2.min(), paras[l].u2.max(), paras[l].v2.min(), paras[l].v2.max());

            warpImage(tmp, paras[l].im1, paras[l].im2, paras[l].u1, paras[l].v1);
            sprintf(buf, out, frameID, iter, l, "sub-warp");
            imwrite(buf, tmp);

            warpImage(tmp, paras[l].im2, paras[l].im1, paras[l].u2, paras[l].v2);
            sprintf(buf, out, frameID, iter, l, "sub-warpr");
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
        segment(centers1, centers2, masks[curIdx], masks[next], numOfLabels,
                im[curIdx], im[next], flow, rflow);
        
        // for test only
        sprintf(buf, out, frameID, iter, 0, "final-layers1");
        flow2color(ucimg, masks[curIdx], masks[curIdx]);
        imwrite(buf, ucimg);

        coverLabels(tmp, im[curIdx], ucimg);
        sprintf(buf, out, frameID, iter, 0, "merge-final-layers1");
        imwrite(buf, tmp);

        sprintf(buf, out, frameID, iter, 0, "final-layers2");
        flow2color(ucimg, masks[next], masks[next]);
        imwrite(buf, ucimg);

        coverLabels(tmp, im[next], ucimg);
        sprintf(buf, out, frameID, iter, 0, "merge-final-layers2");
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

void LAOF::transferLabels(DImage &dstLabels, const DImage &srcLabels, const DImage &flow)
{
    // printf("transfer label... ");

    int width = srcLabels.nWidth(), height = srcLabels.nHeight();
    int channels = srcLabels.nChannels();
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
                dstLabels[noffset+k] = srcLabels[offset+k];
        }
    }

    // printf("done\n");
}

void LAOF::segment(DImage &centers1, DImage &centers2,
                   DImage &layers1, DImage &layers2, int &numOfLabels,
                  const DImage &im1, const DImage &im2,
                  const DImage &flow1, const DImage &flow2)
{
    assert(im1.match3D(im2) && flow1.match3D(flow2) && flow1.nChannels() == 2 &&
           im1.match2D(flow1));

    const int maxNumOfSegs = 10;   
    int width = im1.nWidth(), height = im1.nHeight();
    DImage warpF1, warpF2, warpI1, warpI2, fInfo1, fInfo2, rfInfo, tmp;
    UCImage ucimg;
    MotionLayers ml;
    std::vector<DImage> vec;
    
    // warp flow 2 to flow 1
    warpImage(warpF2, flow2, flow2, flow1);

    ml.flowInfo(fInfo1, flow1);
    ml.flowInfo(rfInfo, warpF2);
    vec.push_back(fInfo1);
    vec.push_back(rfInfo);
    mergew(fInfo1, vec);
        
    // kmeans
    if (centers1.isEmpty())
    {
        numOfLabels = ml.cluster(centers1, layers1, fInfo1, width, height,
                                 2, maxNumOfSegs, 13, true);

        // for test only
        sprintf(buf, out, frameID, 0, 0, "layers1");
        flow2color(ucimg, layers1, layers1);
        imwrite(buf, ucimg);

        coverLabels(tmp, im1, ucimg);
        sprintf(buf, out, frameID, 0, 0, "merge-layers1");
        imwrite(buf, tmp);
    }
    
    // refine layers of im1
    printf("refining layers of im1...\n");
    warpImage(warpI2, im1, im2, flow1);
    ml.refine(centers1, layers1, numOfLabels, im1, warpI2, fInfo1);

    // for test only
    sprintf(buf, out, frameID, 0, 0, "layers1-refine");
    flow2color(ucimg, layers1, layers1);
    imwrite(buf, ucimg);

    coverLabels(tmp, im1, ucimg);
    sprintf(buf, out, frameID, 0, 0, "merge-layers1-refine");
    imwrite(buf, tmp);
    
    // warp flow 2 to flow 1
    printf("cacluating layers of im2...\n");
    warpImage(warpF1, flow1, flow1, flow2);
    ml.flowInfo(fInfo2, warpF1);
    ml.flowInfo(rfInfo, flow2);
    vec.clear();
    vec.push_back(fInfo2);
    vec.push_back(rfInfo);
    mergew(fInfo2, vec);

    // transfer labels to im2
    transferLabels(layers2, layers1, flow1);
    ml.createCenterByLabels(centers2, numOfLabels, layers2, fInfo2);

    // refine layers of im2
    printf("refining layers of im2...\n");
    warpImage(warpI1, im2, im1, flow2);
    ml.refine(centers2, layers2, numOfLabels, warpI1, im2, fInfo2);

    // transfer labels to im1
    transferLabels(layers1, layers2, flow2);
    ml.createCenterByLabels(centers1, numOfLabels, layers1, fInfo1);
    ml.refine(centers1, layers1, numOfLabels, im1, warpI2, fInfo1);
}

void LAOF::intFlow(DImage &dst, const DImage &src)
{
    dst.create(src.nWidth(), src.nHeight(), src.nChannels());

    for (int i = 0; i < src.nElements(); ++i)
    {
        if (src[i] < 0)
            dst[i] = (int)(src[i] - 0.5);
        else
            dst[i] = (int)(src[i] + 0.5);
    }
}
