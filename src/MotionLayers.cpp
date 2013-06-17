#include "GCoptimization.h"
#include "SuperPixels.h"
#include "MotionLayers.h"
#include "ML.h"

#ifdef DEBUG
#include <highgui.h>
using namespace cv;

static int frame = 0;
char buf[256];
const char *pattern = "%s%03d.jpg";
#endif

using cv::Point2d;

const int imHistSize = 12;
const int oriHistSize = 8;
const int magHistSize = 11;

double smoothFn(int p1, int p2, int l1, int l2, void *data)
{
    if (l1 == l2) return 0;

    double *pdata = (double *)data;
    const double sigma2 = 4;
    double pDist = dist2(pdata+p1, pdata+p2, 0, 2);
    double cDist = dist2(pdata+p1, pdata+p2, 2, 5);

    return exp(-(cDist*cDist) / (pDist*2*sigma2));
}

void MotionLayers::init(const Mat &im, const Mat &flow)
{
    assert(im.channels() == 3 && flow.channels() == 2);

    if (im.depth() != CV_IMAGE)
        im.convertTo(m_im, CV_IMAGE);
    else
        im.copyTo(m_im);

    if (flow.depth() != CV_PERCISION)
        flow.convertTo(m_flow, CV_PERCISION);
    else
        flow.copyTo(m_flow);
    
    m_width = im.cols, m_height = im.rows;
}

void MotionLayers::covariance(vector<PERCISION> &cov, vector< vector<Point> > &pos)
{
    int nlabel = pos.size();
    
    // cov(X, Y) = E(XY) - E(X)E(Y)
    cov.resize(nlabel);
    memset(cov.data(), 0, sizeof(PERCISION) * nlabel);

    double EX = 0, EY = 0, EXY = 0;
    for (int l = 0; l < nlabel; ++l)
    {
        int size = pos[l].size();
        for (int i = 0; i < size; ++i)
        {
            Point p = pos[l][i];
            MOTION motion = m_flow.at<MOTION>(p.y, p.x);

            EX += motion[0];
            EY += motion[1];
            EXY += motion[0] * motion[1];
        }

        EX /= size, EY /= size, EXY /= size;
        cov[l] = EXY - EX * EY;
    }
}

// use array of points to represent labels
void MotionLayers::extract(vector< vector<Point> > &pos, int nlabel, const Mat &label)
{
    assert(label.cols == m_width && label.rows == m_height &&
           label.depth() == CV_LABEL);

    pos.resize(nlabel);
    for (int l = 0; l < nlabel; ++l)
        pos[l].clear();

    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int l = label.at<LABEL>(h, w);
            pos[l].push_back(Point(w, h));
        }
}

// patch = spatial info + color histogram + motion histogram
void MotionLayers::extractAsHist(vector<SPPATCH> &patch, const vector< vector<Point> > &pos)
{
    int nlabel = pos.size();
    patch.resize(nlabel);
    
    Mat imPatch, oriPatch, magPatch, lab, hist;
    Point2d point;
    const int channels = 3;
    int channel0[] = {0}, channel1[] = {1}, channel2[] = {2};
    float lRange[] = {0, 101}, aRange[] = {-127, 128}, bRange[] = {-127, 128};
    const float *lRanges[] = {lRange}, *aRanges[] = {aRange}, *bRanges[] = {bRange};

    float oriRange[] = {0, 2 * CV_PI + 1};
    const float *oriRanges[] = {oriRange};
    
    // 假设 u, v的常规范围为 [-10, 10]
    float magRange[magHistSize + 1] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, DBL_MAX};
    const float *magRanges[] = {magRange};

    int ndata = 2 + imHistSize * channels + oriHistSize + magHistSize;
    cvtColor(m_im, lab, CV_BGR2Lab);
    for (int l = 0; l < nlabel; ++l)
    {
        patch[l].resize(ndata);
        
        int size = pos[l].size();
        imPatch.create(1, size, lab.type());
        oriPatch.create(1, size, CV_PERCISION);
        magPatch.create(1, size, CV_PERCISION);
        point.x = 0, point.y = 0;
        for (int i = 0; i < size; ++i)
        {
            Point p = pos[l][i];
            imPatch.at<PIXEL>(i) = lab.at<PIXEL>(p.y, p.x);
            MOTION m = m_flow.at<MOTION>(p.y, p.x);
            oriPatch.at<PERCISION>(i) = atan2(m[1], m[0]) + CV_PI; // [0, 2*pi]
            magPatch.at<PERCISION>(i) = m[1] * m[1] + m[0] * m[0];
            point.x += p.x, point.y += p.y;
        }
        
        // 位置信息
        int idx = 0;
        patch[l][idx++] = point.x / size;
        patch[l][idx++] = point.y / size;
        assert(patch[l][0] < m_width && patch[l][0] >= 0);
        assert(patch[l][1] < m_height && patch[l][1] >= 0);
                                                                  
        // Lab直方图
        calcHist(&lab, 1, channel0, Mat(), hist, 1, &imHistSize, lRanges);
        for (int i = 0; i < imHistSize; ++i, ++idx)
        {
            patch[l][idx] = hist.at<float>(i) / size;
            assert(patch[l][idx] <= 1 && patch[l][idx] >= 0);
        }
        
        calcHist(&lab, 1, channel1, Mat(), hist, 1, &imHistSize, aRanges);
        for (int i = 0; i < imHistSize; ++i, ++idx)
        {
            patch[l][idx] = hist.at<float>(i) / size;
            assert(patch[l][idx] <= 1 && patch[l][idx] >= 0);
        }

        calcHist(&lab, 1, channel2, Mat(), hist, 1, &imHistSize, bRanges);
        for (int i = 0; i < imHistSize; ++i, ++idx)
        {
            patch[l][idx] = hist.at<float>(i) / size;
            assert(patch[l][idx] <= 1 && patch[l][idx] >= 0);
        }
        
        // 运动方向直方图
        calcHist(&oriPatch, 1, channel0, Mat(), hist, 1, &oriHistSize, oriRanges);
        for (int i = 0; i < oriHistSize; ++i, ++idx)
        {
            patch[l][idx] = hist.at<float>(i) / size;
            assert(patch[l][idx] <= 1 && patch[l][idx] >= 0);
        }

        // 运动幅度直方图
        calcHist(&magPatch, 1, channel0, Mat(), hist, 1, &magHistSize, magRanges, false);
        for (int i = 0; i < magHistSize; ++i, ++idx)
        {
            patch[l][idx] = hist.at<float>(i) / size;
            assert(patch[l][idx] <= 1 && patch[l][idx] >= 0);
        }
    }
}

void MotionLayers::extractAsPixel(vector<SPPATCH> &patch, const vector< vector<Point> > &pos)
{
    int ndata = 2 + 3 + 2, nlabel = pos.size();
    patch.resize(nlabel);

    for (int i = 0; i < nlabel; ++i)
    {
        patch[i].resize(ndata);
        memset(patch[i].data(), 0, sizeof(PERCISION) * ndata);
    }

    for (int l = 0; l < nlabel; ++l)
    {
        int size = pos[l].size();
        for (int i = 0; i < size; ++i)
        {
            Point p = pos[l][i];
            PIXEL pixel = m_im.at<PIXEL>(p.y, p.x);
            MOTION motion = m_flow.at<MOTION>(p.y, p.x);

            patch[l][0] += p.x;
            patch[l][1] += p.y;
            patch[l][2] += pixel[0];
            patch[l][3] += pixel[1];
            patch[l][4] += pixel[2];
            patch[l][5] += motion[0];
            patch[l][6] += motion[1];
        }

        for (int i = 0; i < ndata; ++i)
            patch[l][i] /= size;
    }
}

void MotionLayers::initSegment(int ncluster)
{
    // 计算 super pixel
    printf("generate superpixels...\n");
    SuperPixels sp(50);
    sp.setSourceImage(m_im);
    sp.generate();
    m_nSPLbl = sp.getLabels(m_spLbl);
    printf("done, get %d super pixels\n", m_nSPLbl);

#ifdef DEBUG
    sprintf(buf, pattern, "superpixels", frame);
    sp.saveSegmentationEdges(buf);
#endif

    // 用直方图表示 super pixel
    printf("represent super pixels...\n");
    extract(m_sp2pos, m_nSPLbl, m_spLbl);
    extractAsPixel(m_spPixel, m_sp2pos);
    extractAsHist(m_spHist, m_sp2pos);
    printf("done\n");

#ifdef DEBUG
    Mat img(m_height, m_width, m_im.type());
    PIXEL pixel;
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int lbl = m_spLbl.at<LABEL>(h, w);
            pixel[0] = m_spPixel[lbl][2];
            pixel[1] = m_spPixel[lbl][3];
            pixel[2] = m_spPixel[lbl][4];
            img.at<PIXEL>(h, w) = pixel;
        }

    sprintf(buf, pattern, "spAsColor", frame);
    imwrite(buf, img);

    double maxu, maxv, minu, minv, maxrad, u, v;
    maxu = minu = m_spPixel[0][5];
    maxv = minv = m_spPixel[0][6];
    maxrad = minu*minu + minv*minv;
    for (int i = 1; i < m_nSPLbl; ++i)
    {
        u = m_spPixel[i][5], v = m_spPixel[i][6];
        maxu = std::max(maxu, u);
        minu = std::min(minu, u);
        maxv = std::max(maxv, v);
        minv = std::min(minv, v);
        maxrad = std::max(maxrad, u*u+v*v);
    }
    printf("u: %.6f .. %.6f, v: %.6f .. %.6f\n", minu, maxu, minv, maxv);
    printf("maxrad: %.6f\n", maxrad);
    
    double ufactor, vfactor, radfactor;
    ufactor = 255. / (maxu - minu);
    vfactor = 255. / (maxv - minv);
    radfactor = 255. / maxrad;
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int lbl = m_spLbl.at<LABEL>(h, w);
            u = m_spPixel[lbl][5];
            v = m_spPixel[lbl][6];
            pixel[0] = (u - minu) * ufactor;
            pixel[1] = (v - minv) * vfactor;
            pixel[2] = (u*u + v*v) * radfactor;
            img.at<PIXEL>(h, w) = pixel;
        }

    sprintf(buf, pattern, "spAsMotion", frame);
    imwrite(buf, img);

    int bin_w = 10, bin_h = 100;
    int histBins = m_spHist[0].size();
    int color = 255. / (histBins+1);
    for (int l = 0; l < m_nSPLbl; ++l)
    {
        Mat histImg = Mat::zeros(bin_h, bin_w * histBins, CV_8UC3);
        for(int i = 2; i < histBins; ++i)
        {
            PERCISION val = m_spHist[l][i];
            printf("val: %.6f\n", val);
            assert(val >= 0 && val <= 1);
            
            Point lbp(i * bin_w, 0);
            Point rup((i+1) * bin_w, cvRound(val * bin_h));
            rectangle(histImg, lbp, rup, Scalar(color*(i+1), color*(i-1), color*i), -1);
        }
        sprintf(buf, pattern, "spAsHist", frame*100+l);
        imwrite(buf, histImg);
    }
    
#endif
    
    // 计算可信的 patch
    covariance(m_cov, m_sp2pos);

#ifdef DEBUG
    double minCov = DBL_MAX, maxCov = 0;
    for (int i = 0; i < m_nSPLbl; ++i)
    {
        minCov = std::min(minCov, fabs(m_cov[i]));
        maxCov = std::max(maxCov, fabs(m_cov[i]));
    }

    printf("cov: %.6f .. %.6f\n", minCov, maxCov);
#endif
    
    const PERCISION threshold = 5e-3;
    vector<SPPATCH> samples;
    vector<int> spid2tid(m_nSPLbl, -1);
    int trust = 0;
    for (int l = 0; l < m_nSPLbl; ++l)
    {
        if (fabs(m_cov[l]) < threshold)
        {
            samples.push_back(m_spPixel[l]);
            spid2tid[l] = trust++;
        }
    }
    printf("trusted points: %d / %d\n", trust, m_nSPLbl);
    
#ifdef DEBUG
    Mat trustImg = Mat::zeros(m_height, m_width, m_im.type());
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int lbl = m_spLbl.at<LABEL>(h, w);
            
            if (spid2tid[lbl] >= 0)
            {
                PIXEL color;
                color[0] = m_spPixel[lbl][2];
                color[1] = m_spPixel[lbl][3];
                color[2] = m_spPixel[lbl][4];

                trustImg.at<PIXEL>(h, w) = color;
            }
        }

    sprintf(buf, pattern, "trust", frame);
    imwrite(buf, trustImg);
#endif
    
    // kmeans
    vector<SPPATCH> centers;
    vector<LABEL> label;
    kmeans(centers, label, samples, ncluster, kmdist);
    printf("kmeans done\n");

    m_label = Mat::zeros(m_height, m_width, CV_LABEL);
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int spLbl = m_spLbl.at<LABEL>(h, w);
            int nLbl = spid2tid[spLbl];
            if (nLbl >= 0)
                m_label.at<LABEL>(h, w) = label[nLbl] + 1;
        }

#ifdef DEBUG
    sprintf(buf, pattern, "kmeans", frame);
    imwrite(buf, m_label * 255 / (ncluster+1));
#endif
}

void MotionLayers::refineSegment(int ncluster)
{
    // prepare data for graph cut
    printf("prepare data for graph cut...\n");
    int nPixelLen = m_spPixel[0].size();
    vector<double> spPixel(m_nSPLbl*nPixelLen);

    for (int i = 0; i < m_nSPLbl; ++i)
        for (int l = 0; l < nPixelLen; ++l)
            spPixel[i*nPixelLen+l] = m_spPixel[i][l];

    printf("done\n");
    
    // extract init layers
    printf("extract init layers\n");
    vector< vector<Point> > npos;    
    vector<SPPATCH> layer;
    extract(npos, ncluster+1, m_label);
    extractAsHist(layer, npos);
    printf("done\n");
    
    // prepare data term
    printf("prepare data term\n");
    int nHistLen = layer[0].size();
    vector<double> data(ncluster*m_nSPLbl);

    PERCISION weight, kl;
    for (int i = 0; i < m_nSPLbl; ++i)
    {
        weight = 1. / sqrt(1 + m_cov[i] * m_cov[i]);
        for (int l = 1; l <= ncluster; ++l)
        {
            kl = skldist(m_spHist[i].data(), layer[l].data(), 0, nHistLen);
            data[i*ncluster+l-1] = weight * kl;
        }
    }
    printf("done\n");

    printf("running graph cut, sites: %d, labels: %d...\n", m_nSPLbl, ncluster);
    GCoptimizationGeneralGraph gc(m_nSPLbl, ncluster);

    printf("set data cost...\n");
    gc.setDataCost(data.data());
    gc.setSmoothCost(smoothFn, spPixel.data());
    printf("done\n");

    // set up graph
    // horizontal
    printf("set up horizontal connection...\n");
    for (int h = 0; h < m_height; ++h)
        for (int w = 1; w < m_width; ++w)
        {
            int lbl = m_spLbl.at<LABEL>(h, w);
            int lblx_1 = m_spLbl.at<LABEL>(h, w-1);
            assert(lbl >= 0 && lblx_1 >= 0 && m_nSPLbl > lbl && m_nSPLbl > lblx_1);

            if (lbl != lblx_1)
                gc.setNeighbors(lbl, lblx_1);
        }

    // vertical
    printf("set up vertical connection...\n");
    for (int h = 1; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int lbl = m_spLbl.at<LABEL>(h, w);
            int lbly_1 = m_spLbl.at<LABEL>(h-1, w);
            assert(lbl >= 0 && lbly_1 >= 0 && m_nSPLbl > lbl && m_nSPLbl > lbly_1);
                
            if (lbl != lbly_1)
                gc.setNeighbors(lbl, lbly_1);
        }
    printf("done\n");
        
    printf("Before optimization energy is %.6f\n",gc.compute_energy());
    gc.expansion(2);
    printf("After optimization energy is %.6f\n",gc.compute_energy());

    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int spLbl = m_spLbl.at<LABEL>(h, w);
            m_label.at<LABEL>(h, w) = gc.whatLabel(spLbl);
        }

#ifdef DEBUG
    sprintf(buf, pattern, "refine", frame++);
    imwrite(buf, m_label * 255 / ncluster);
#endif
}

// distance for kmeans
double MotionLayers::kmdist(const PERCISION *p1, const PERCISION *p2, int start, int end)
{
    const double spatial = 0.5;
    const double color = 0;
    const double motion = 1.;
    double pDist, cDist, oDist, mDist;
    int idx = start;
    
    // spatial info
    pDist = dist2(p1, p2, idx, idx+2);
    pDist = log10(pDist+9);
    idx += 2;

/*    // color info
    cDist = skldist(p1, p2, idx, idx+imHistSize*3);
    idx += imHistSize*3;
    
    // flow info
    oDist = skldist(p1, p2, idx, idx+oriHistSize);
    idx += oriHistSize;

    mDist = skldist(p1, p2, idx, idx+magHistSize);
    idx += magHistSize;

    assert(idx <= end);
    return (spatial * pDist + color * cDist + motion * (oDist + mDist));
*/
    cDist = dist2(p1, p2, idx, idx+3);
    idx += 3;

    mDist = dist2(p1, p2, idx, idx+2);
    idx += 2;
    assert(idx <= end);

    return (spatial * pDist + color * cDist + motion * mDist);
}

// kullback-leiber divergence
double MotionLayers::kldist(const PERCISION *p1, const PERCISION *p2, int start, int end)
{
    double cost = 0;
    for (int i = start; i < end; ++i)
        if (p1[i] > ESP && p2[i] > ESP)
            cost += log(p1[i] / p2[i]) * p1[i];

    return cost;
}

// symmetric kullback-leiber
double MotionLayers::skldist(const PERCISION *p1, const PERCISION *p2, int start, int end)
{
    return kldist(p1, p2, start, end) + kldist(p2, p1, start, end);
}
