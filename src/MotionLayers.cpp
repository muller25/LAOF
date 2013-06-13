#include "GCoptimization.h"
#include "SuperPixels.h"
#include "MotionLayers.h"
#include "ML.h"

using cv::Point2d;

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

    if (im.type() != CV_PERCSION)
        im.convertTo(m_im, CV_PERCSION);
    else
        im.copyTo(m_im);

    if (flow.type() != CV_PERCSION)
        flow.convertTo(m_flow, CV_PERCSION);
    else
        flow.copyTo(m_flow);
    
    m_width = im.cols, m_height = im.rows;
}

void MotionLayers::covariance(vector<PERCISION> &cov, vector< vector<Point> > &pos)
{
    // cov(X, Y) = E(XY) - E(X)E(Y)
    cov.resize(m_nSPLbl);
    memset(cov.data(), 0, sizeof(PERCISION) * m_nSPLbl);

    PERCISION EX = 0, EY = 0, EXY = 0;
    for (int l = 0; l < m_nSPLbl; ++l)
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

void MotionLayers::extract(vector< vector<Point> > &pos, int nlabel, const Mat &label)
{
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
    const int channels = 3, channel = 1;
    float lRange[] = {0, 101}, aRange[] = {-127, 128}, bRange[] = {-127, 128};
    const float *lRanges[] = {lRange}, *aRanges[] = {aRange}, *bRanges[] = {bRange};

    float oriRange[] = {0, 361};
    const float *oriRanges[] = {oriRange};
    
    // 假设 u, v的常规范围为 [-10, 10]
    float magRange[magHistSize + 1] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, DBL_MAX};
    const float *magRanges[] = {magRange};
    
    cvtColor(m_im, lab, CV_BGR2Lab);
    for (int l = 0; l < nlabel; ++l)
    {
        patch[l].resize(2 + imHistSize * channels + oriHistSize + magHistSize);
        
        int size = pos[l].size();
        imPatch.create(1, size, m_im.depth());
        oriPatch.create(1, size, CV_PERCSION);
        magPatch.create(1, size, CV_PERCSION);
        point.x = 0, point.y = 0;
        for (int i = 0; i < size; ++i)
        {
            Point p = pos[l][i];
            imPatch.at<PIXEL>(i) = lab.at<PIXEL>(p.y, p.x);
            MOTION m = m_flow.at<MOTION>(p.y, p.x);
            oriPatch.at<PERCISION>(i) = atan2(m[1], m[0]) + CV_PI; // [0, 2*pi]
            magPatch.at<PERCISION>(i) = m[1] * m[1] + m[0] * m[0];
            point.x += p.x;
            point.y += p.y;
        }

        // 位置信息
        int idx = 0;
        patch[l][idx++] = point.x / size;
        patch[l][idx++] = point.y / size;
        
        // Lab直方图
        int c = 0;
        calcHist(&lab, 1, &c, Mat(), hist, 1, &imHistSize, lRanges);
        for (int i = 0; i < imHistSize; ++i)
            patch[l][idx++] = hist.at<float>(i) / size;

        c = 1;
        calcHist(&lab, 1, &c, Mat(), hist, 1, &imHistSize, aRanges);
        for (int i = 0; i < imHistSize; ++i)
            patch[l][idx++] = hist.at<float>(i) / size;

        c = 2;
        calcHist(&lab, 1, &c, Mat(), hist, 1, &imHistSize, bRanges);
        for (int i = 0; i < imHistSize; ++i)
            patch[l][idx++] = hist.at<float>(i) / size;
        
        // 运动方向直方图
        calcHist(&oriPatch, 1, &channel, Mat(), hist, 1, &oriHistSize, oriRanges);
        for (int i = 0; i < oriHistSize; ++i)
            patch[l][idx++] = hist.at<float>(i) / size;

        // 运动幅度直方图
        calcHist(&magPatch, 1, &channel, Mat(), hist, 1, &magHistSize, magRanges, false);
        for (int i = 0; i < magHistSize; ++i)
            patch[l][idx++] = hist.at<float>(i) / size;
    }
}

void MotionLayers::extractAsPixel(vector<SPPATCH> &patch, const vector< vector<Point> > &pos)
{
    int ndata = 2 + 3 + 2;
    int nlabel = pos.size();
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
    SuperPixels sp;
    sp.setSourceImage(m_im);
    sp.generate();
    m_nSPLbl = sp.getLabels(m_spLbl);

    // 用直方图表示 super pixel
    extract(m_sp2pos, m_nSPLbl, m_spLbl);
    extractAsHist(m_spHist, m_sp2pos);

    // 计算可信的 patch
    covariance(m_cov, m_sp2pos);

    vector<SPPATCH> samples, centers;
    vector<LABEL> label;
    vector<int> spid2tid(m_nSPLbl, -1);
    int trust = 0;
    for (int l = 0; l < m_nSPLbl; ++l)
        if (abs(m_cov[l]) < ESP)
        {
            samples.push_back(m_spHist[l]);
            spid2tid[l] = trust++;
        }

    // kmeans
    kmeans(centers, label, samples, ncluster, kmdist);
    m_label = Mat::zeros(m_height, m_width, CV_LABEL);
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int spLbl = m_spLbl.at<LABEL>(h, w);
            int nLbl = spid2tid[spLbl];
            if (nLbl >= 0)
                m_label.at<LABEL>(h, w) = label[nLbl] + 1;
        }
}

void MotionLayers::refineSegment(int ncluster)
{
    vector<SPPATCH> sp, layer;
    vector< vector<Point> > npos;
    
    extractAsPixel(sp, m_sp2pos);
    int nPixelLen = sp[0].size();
    double *spPixel = new double[m_nSPLbl*nPixelLen];
    for (int i = 0; i < m_nSPLbl; ++i)
        for (int l = 0; l < nPixelLen; ++l)
            spPixel[i*nPixelLen+l] = sp[i][l];
    
    extract(npos, ncluster+1, m_label);
    extractAsHist(layer, npos);
    
    int nHistLen = m_spHist[0].size();
    double *data = new double[ncluster * m_nSPLbl];
    PERCISION weight, kl;
    for (int i = 0; i < m_nSPLbl; ++i)
        for (int l = 1; l <= ncluster; ++l)
        {
            weight = 1. / sqrt(1 + m_cov[i] * m_cov[i]);
            kl = skldist(m_spHist[i].data(), layer[l].data(), 0, nHistLen);
            data[i*m_nSPLbl+l-1] = weight * kl;
        }

    try {
        GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(m_nSPLbl, ncluster);
		gc->setDataCost(data);
        gc->setSmoothCost(smoothFn, spPixel);

        // horizontal
        for (int h = 0; h < m_height; ++h)
            for (int w = 1; w < m_width; ++w)
            {
                int lbl = m_spLbl.at<LABEL>(h, w);
                int lblx_1 = m_spLbl.at<LABEL>(h, w-1);
                if (lbl != lblx_1)
                    gc->setNeighbors(lbl, lblx_1);
            }

        // vertical
        for (int h = 1; h < m_height; ++h)
            for (int w = 0; w < m_width; ++w)
            {
                int lbl = m_spLbl.at<LABEL>(h, w);
                int lbly_1 = m_spLbl.at<LABEL>(h-1, w);
                if (lbl != lbly_1)
                    gc->setNeighbors(lbl, lbly_1);
            }
        
		printf("Before optimization energy is %.6f\n",gc->compute_energy());
		gc->expansion(2);
		printf("After optimization energy is %.6f\n",gc->compute_energy());

        for (int h = 0; h < m_height; ++h)
            for (int w = 0; w < m_width; ++w)
            {
                int spLbl = m_spLbl.at<LABEL>(h, w);
                m_label.at<LABEL>(h, w) = gc->whatLabel(spLbl);
            }

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}
    
    delete []data;
    delete []spPixel;
}

// distance for kmeans
double MotionLayers::kmdist(const PERCISION *p1, const PERCISION *p2, int start, int end)
{
    double spatial = 1.;
    double color = 0;
    double motion = 1.;
    double pDist, cDist, oDist, mDist;
    int idx = start;
    
    // spatial info
    pDist = dist2(p1, p2, idx, idx+2);
    pDist = log10(pDist + 9);
    idx += 2;

    // color info
    cDist = kldist(p1, p2, idx, idx+imHistSize*3) + kldist(p2, p1, idx, idx+imHistSize*3);
    idx += imHistSize*3;
    
    // flow info
    oDist = kldist(p1, p2, idx, idx+oriHistSize) + kldist(p2, p1, idx, idx+oriHistSize);
    idx += oriHistSize;

    mDist = kldist(p1, p2, idx, idx+magHistSize) + kldist(p2, p1, idx, idx+magHistSize);
    idx += magHistSize;
    
    return (spatial * pDist + color * cDist + motion * (oDist + mDist));
}

// kullback-leiber divergence
double MotionLayers::kldist(const PERCISION *p1, const PERCISION *p2, int start, int end)
{
    PERCISION cost = 0;

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
