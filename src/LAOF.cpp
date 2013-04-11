#include "LAOF.h"
#include "OpticalFlow.h"
#include "MotionLayers.h"
#include "ImageIO.h"
#include "Flow2Color.h"
#include "OFPara.h"
#include "Utils.h"

// init optical flow parameters
const double as = 0.026;
const double ap = 0.012;
const double ratio = 0.75;
const int minWidth = 20;
const int nBiIter = 5;//14;
const int nIRLSIter = 1;
const int nSORIter = 10;//30;

extern int frameStart;
extern int nsegs;

char buf[256];
char LAOF::out[256];

void LAOF::EM(std::vector<DImage> &u, std::vector<DImage> &v,
              std::vector<DImage> &ur, std::vector<DImage> &vr,
              std::vector<DImage> &mask, const std::vector<DImage> &im,
              int cur, int nIters)
{
    assert(u.size() == v.size() && v.size() == ur.size() &&
           ur.size() == mask.size() && mask.size() == im.size() && im.size() == 2);

    int next = 1 - cur;
    int width = im[0].nWidth(), height = im[0].nHeight(), nlabels = 1;
    DImage tmp, centers1, centers2;
    UCImage ucimg;
    std::vector<OFPara> paras;
    OpticalFlow of;
    
    // init mask
    mask[cur].create(width, height);
    mask[next].create(width, height);
    
    // start EM iterations
    for (int iter = 0; iter < nIters; ++iter)
    {
        printf("EM iterations: %d\n", iter);

        // split images, mask into pieces
        OFPara::split(paras, nlabels, im[cur], im[next], mask[cur], mask[next]);

        // loop over layers and run biDir optical flow
        printf("loop over %d layers and run biDir optical flow\n", nlabels);
        for (int l = 0; l < nlabels; ++l)
        {
            // for test only
            if (l > 0)
            {
                sprintf(buf, out, frameStart, iter, l, "mask1");
                imwrite(buf, paras[l].mask1);
            
                sprintf(buf, out, frameStart, iter, l, "mask2");
                imwrite(buf, paras[l].mask2);
            }
            
            of.biC2FFlow(paras[l].u1, paras[l].v1, paras[l].u2, paras[l].v2,
                         paras[l].im1, paras[l].im2,
                         paras[l].mask1, paras[l].mask2,
                         as, ap, ratio, minWidth, nBiIter, nIRLSIter, nSORIter);

            // for test only
            printf("u : %.6f .. %.6f, v : %.6f .. %.6f\n", paras[l].u1.min(), paras[l].u1.max(), paras[l].v1.min(), paras[l].v1.max());
            printf("ur: %.6f .. %.6f, vr: %.6f .. %.6f\n", paras[l].u2.min(), paras[l].u2.max(), paras[l].v2.min(), paras[l].v2.max());
        }

        OFPara::restore(u[cur], v[cur], ur[cur], vr[cur], paras, width, height);

        // for test only
        printf("****** global flow ******\n");
        sprintf(buf, out, frameStart, iter, 0, "flow");
        flow2color(ucimg, u[cur], v[cur]);
        imwrite(buf, ucimg);

        sprintf(buf, out, frameStart, iter, 0, "rflow");
        flow2color(ucimg, ur[cur], vr[cur]);
        imwrite(buf, ucimg);

        // extract motion layers
        printf("****** extract motion layer ******\n");
        nlabels = nsegs;
        segment(centers1, centers2, mask[cur], mask[next], nlabels,
                im[cur], im[next], u[cur], v[cur], ur[cur], vr[cur]);
    }

    frameStart++;
}

void LAOF::segment(DImage &centers1, DImage &centers2,
                   DImage &layers1, DImage &layers2, int &nlabels,
                   const DImage &im1, const DImage &im2,
                   const DImage &u, const DImage &v,
                   const DImage &ur, const DImage &vr)
{
    DImage warpI2, warpI1, tmp;
    UCImage ucimg;
    MotionLayers ml;

    warpImage(warpI2, im1, im2, u, v);
    warpImage(warpI1, im2, im1, ur, vr);
    
    // kmeans
    ml.kcluster(centers1, layers1, nlabels, im1, warpI2, u, v);

    // for test only
    showLayers(ucimg, layers1);
    sprintf(buf, out, frameStart, 0, 0, "layers1");
    imwrite(buf, ucimg);

    // show image under each layer
    for (int l = 0; l < nlabels; ++l)
    {
        genLayerMask(ucimg, layers1, l);
        cut(tmp, im1, ucimg);
        sprintf(buf, out, frameStart, 0, l, "kseg");
        imwrite(buf, tmp);
    }
    
    // refine layers of im1
    printf("refining layers of im1...\n");
    ml.refine(centers1, layers1, nlabels, im1, warpI2, u, v);

    
    ml.kcluster(centers2, layers2, nlabels, warpI1, im2, ur, vr);
    ml.refine(centers2, layers2, nlabels, warpI1, im2, ur, vr);

    // for test only
    showLayers(ucimg, layers2);
    sprintf(buf, out, frameStart, 0, 0, "layers2");
    imwrite(buf, ucimg);

    // show image under each layer
    for (int l = 0; l < nlabels; ++l)
    {
        genLayerMask(ucimg, layers2, l);
        cut(tmp, im1, ucimg);
        sprintf(buf, out, frameStart, 1, l, "kseg");
        imwrite(buf, tmp);
    }
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
