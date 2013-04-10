#ifndef _MotionLayers_H
#define _MotionLayers_H

#include "Image.h"
#include "ImageIO.h"
#include <highgui.h>

class MotionLayers
{
public:
    void spatialInfo(DImage &info, const DImage &im);
    void imInfo(DImage &info, const DImage &im);
    void flowInfo(DImage &info, const DImage &flow);

    int cluster(DImage &centers, DImage &layers, const DImage &im,
                const DImage &flow, const DImage &rflow,
                int start=2, int end=10, double na=15, bool reArrange=false);

    int cluster(DImage &centers, DImage &layers,
                const DImage &features, int width, int height,
                int start=2, int end=10, double na=15, bool reArrange=false);

    void refine(DImage &centers, DImage &layers, int labels,
                const DImage &im1, const DImage &im2,
                const DImage &features);

    void kcluster(DImage &centers, DImage &layers, int nlabels,
                  const DImage &u, const DImage &v);

    void scluster(DImage &centers, DImage &layers, int numOfClusters,
                  const DImage &u, const DImage &v);

    void refine(DImage &centers, DImage &layers, int numOfClusters,
                const DImage &im1, const DImage &im2,
                const DImage &u, const DImage &v);

    template <class I, class M>
    void LabComfirmity(DImage &prob, const Image<I> &im, const Image<M> &mask);

    template <class M>
    void OMComfirmity(DImage &prob, const DImage &om, const Image<M> &mask);

    // make sure layer order
    template <class T>
    void reArrangeLabels(Image<T> &layers, int labels);

    // create cluster centers by labels
    template <class T, class T1>
    void createCenterByLabels(Image<T> &centers, int numOfLabels,
                              const Image<T1> &labels, const Image<T> &samples);

    inline void dataFn(double *data, int nlabels, const DImage &layers,
                       const DImage &u, const DImage &v, const DImage &im);
    static double smoothFn(int p1, int p2, int l1, int l2, void *pData);
    inline static double kmdist(double *p1, double *p2, int start, int end);

private:
    const static int sWidth = 2;// spatial info width
    const static int iWidth = 3;// image info width
    const static int fWidth = 2;// flow info width
};

// arrange layers from largest component to smallest
template <class T>
void MotionLayers::reArrangeLabels(Image<T> &layers, int labels)
{
    assert(layers.nChannels() == 1);
    
    int *count = new int[labels];
    int *mapping = new int[labels];
    int idx, maxNum, maxid, id;

    memset(count, 0, sizeof(int) * labels);
    for (int i = 0; i < layers.nElements(); ++i)
    {
        idx = layers[i] + 0.5;
        count[idx]++;
    }

    id = 0, maxid = -1;
    for (int i = 0; i < labels; ++i)
    {
        maxNum = -1;
        for (int j = 0; j < labels; ++j)
        {
            if (maxNum < count[j]) {
                maxNum = count[j];
                maxid = j;
            }
        }

        mapping[maxid] = id++;
        count[maxid] = -1;
    }

    for (int i = 0; i < layers.nElements(); ++i)
    {
        idx = layers[i] + 0.5;
        layers[i] = mapping[idx];
    }
    
    delete []mapping;
    delete []count;
}

template <class T, class T1>
void MotionLayers::createCenterByLabels(Image<T> &centers, int numOfLabels,
                                        const Image<T1> &labels, const Image<T> &samples)
{
    // printf("create center by label... ");
    
    int *count = new int[numOfLabels];
    int width = samples.nWidth(), lsize = labels.nSize(), idx;
    
    centers.create(width, numOfLabels);
    memset(count, 0, sizeof(int)*numOfLabels);

    for (int i = 0; i < lsize; ++i)
    {
        idx = labels[i];
        count[idx]++;
        for (int w = 0; w < width; ++w)
            centers[idx*width+w] += samples[i*width+w];
    }

    for (int i = 0; i < numOfLabels; ++i)
        for (int w = 0; w < width; ++w)
            centers[i*width+w] /= (double)count[i];
    
    delete []count;
    // printf("done\n");
}

// 计算mask覆盖的图像区域的lab直方图Hlab
// 然后通过back project的方法计算图像im中每个像素点在Hlab中的概率，得到color comfirmity图像prob
template <class I, class M>
void MotionLayers::LabComfirmity(DImage &prob, const Image<I> &im, const Image<M> &mask)
{
    // convert to lab color space
    cv::Mat imat, lab, maskMat, tmp;

    im.convertTo(tmp);
    tmp.convertTo(imat, CV_32FC3);
    cv::cvtColor(imat, lab, CV_BGR2Lab);

    // convert mask
    mask.convertTo(tmp);
    tmp.convertTo(maskMat, CV_8U);

    const int lbins = 10;
    const int abins = 12;
    const int bbins = 12;
    const int histSize[] = {lbins, abins, bbins};
    const int channels[] = {0, 1, 2};

    // l varies from [0, 100]
    const float lranges[] = {0, 101};
                             
    // a varies from [-127, 127]
    const float aranges[] = {-127, 128};

    // b varies from [-127, 127]
    const float branges[] = {-127, 128};

    const float *ranges[] = {lranges, aranges, branges};
        
    cv::Mat hist;
    cv::calcHist(&imat, 1, channels, maskMat, hist, 3, histSize, ranges, true, false);

    // normalize
    int totalPixels = mask.nonZeros();
    for (int l = 0; l < lbins; ++l)
        for (int a = 0; a < abins; ++a)
            for (int b = 0; b < bbins; ++b)
                hist.at<float>(l, a, b) /= totalPixels;

    // back project
    cv::Mat backproj;
    calcBackProject(&imat, 1, channels, hist, backproj, ranges, 1, true);
    prob.convertFrom(backproj);
    
    // show histogram image
    static int hcount = 0;
    static int ccount = 0;
    static int mcount = 0;
    char buf[256];
    int scale = 10;
    cv::Mat histImgMat = cv::Mat::zeros(abins*scale, bbins*scale, CV_8UC3);
    for (int a = 0; a < abins; ++a)
    {
        for (int b = 0; b < bbins; ++b)
        {
            float c1 = 0, c2 = 0, c3 = 0;
            for (int c = 0; c < lbins; c += 3)
            {
                c1 += hist.at<c, a, b>();
                c2 += hist.at<c+1, a, b>();
                c3 += hist.at<c+2, a, b>();
            }
            
            cv::Scalar color(cvRound(c1*255), cvRound(c2*255), cvRound(c3*255));
            cv::rectangle(histImgMat, cv::Point(a*scale, b*scale),
                          cv::Point((a+1)*scale-1, (b+1)*scale-1),
                          color, CV_FILLED);
        }
    }

    sprintf(buf, "lab-hist-%d.jpg", hcount++);
    cv::imwrite(buf, histImgMat);

    DImage area;
    cut(area, im, mask);
    sprintf(buf, "mask-%d.jpg", mcount++);
    imwrite(buf, area);
    
    sprintf(buf, "lab-comfirmity-%d.jpg", ccount++);
    imwrite(buf, prob);
}

// om为运动u，v的方向-强度图
// 计算mask覆盖的om区域的om直方图Hom
// 然后通过back project的方法计算om中每个元素在Hom中的概率，得到motion comfirmity图像prob
template <class M>
void MotionLayers::OMComfirmity(DImage &prob, const DImage &om, const Image<M> &mask)
{
    // convert om
    cv::Mat omMat, maskMat, tmp;

    om.convertTo(tmp);
    tmp.convertTo(omMat, CV_32FC3);

    // convert mask
    mask.convertTo(tmp);
    tmp.convertTo(maskMat, CV_8U);

    const int obins = 10;
    const int mbins = 10;
    const int histSize[] = {obins, mbins};
    const int channels[] = {0, 1};

    // orientation: [-PI, PI]
    const float oranges[] = {-PI, PI+1};

    // magnitude: [0, 1]
    const float mranges[] = {0, 2};
    const float *ranges[] = {oranges, mranges};
        
    cv::Mat hist;
    cv::calcHist(&omMat, 1, channels, maskMat, hist, 2, histSize, ranges, true, false);

    // probability of each bin
    int totalPixels = mask.nonZeros();
    for (int o = 0; o < obins; ++o)
        for (int m = 0; m < mbins; ++m)
            hist.at<float>(o, m) /= totalPixels;

    // back project
    cv::Mat backproj;
    calcBackProject(&omMat, 1, channels, hist, backproj, ranges, 1, true);
    prob.convertFrom(backproj);
   
    // show image
    static int hcount = 0;
    static int ccount = 0;
    char buf[256];
    int scale = 10;
    cv::Mat histImgMat = cv::Mat::zeros(obins*scale, mbins*scale, CV_8U);
    for (int o = 0; o < obins; ++o)
    {
        for (int m = 0; m < mbins; ++m)
        {
            float binVal = hist.at<float>(o, m);
            int intensity = cvRound(binVal*255);
            cv::rectangle(histImgMat, cv::Point(o*scale, m*scale),
                          cv::Point((o+1)*scale-1, (m+1)*scale-1),
                          intensity, CV_FILLED);
        }
    }

    sprintf(buf, "om-hist-%d.jpg", hcount++);
    cv::imwrite(buf, histImgMat);

    sprintf(buf, "om-comfirmity-%d.jpg", ccount++);
    imwrite(buf, prob);
}

#endif
