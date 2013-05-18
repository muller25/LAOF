#include "RBF.h"

#include <stdio.h>

#include <highgui.h>
using namespace cv;

// #define DEBUG

/*
  描述：
  计算RBF的中心点，中心点必须符合
  1. 梯度幅度大于等于阈值
  2. 相邻中心点的距离大于等于最小距离
  参数：
  centers: 返回值，每个元素由 (x, y) 组成
  gm: 梯度幅度
*/
void RBF::rbf_center(vector<Point> &centers, const Mat &gm)
{
    centers.clear();

    double minVal, maxVal, threshold;
    minMaxLoc(gm, &minVal, &maxVal);
    threshold = maxVal * 0.5;

#ifdef DEBUG
    printf("magnitude: %.6f .. %.6f\n", minVal, maxVal);
    imshow("edges", (gm >= threshold));
    printf("finding centers...\n");
#endif

    int width = gm.cols, height = gm.rows;
    int range = std::max(std::min(width / 10., height / 10.), 2.);
    Mat taken = Mat::zeros(height, width, CV_8U);
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            // 中心点距离测试
            if (taken.at<uchar>(h, w) == 1) continue;
            
            double mag = gm.at<PERCISION>(h, w);
            if (mag < threshold) continue;

            centers.push_back(Point(w, h));
            circle(taken, Point(w, h), range, Scalar(1), -1);
        }
#ifdef DEBUG
    imshow("taken", taken * 255);
    waitKey(0);
#endif
}

void RBF::rbf_solver(Mat &rbfres, const vector<Point> &centers, const Mat &gm)
{
    // estimate weights
    int len = centers.size();
    Mat PHI(len, len, CV_PERCISION), Y(len, 1, CV_PERCISION);

    for (int h = 0; h < len; ++h)
    {
        int x = centers[h].x;
        int y = centers[h].y;
        Y.at<PERCISION>(h, 0) = gm.at<PERCISION>(y, x);
        for (int w = 0; w < len; ++w)
        {
            PERCISION dx = centers[w].x - x;
            PERCISION dy = centers[w].y - y;
            PERCISION dist = sqrt(dx * dx + dy * dy);
            PHI.at<PERCISION>(h, w) = phi(dist);
        }
    }

    Mat PHIt = PHI.t();
    Mat tmp = PHIt * PHI;
    Mat W =  tmp.inv(DECOMP_SVD) * PHIt * Y;

    // interpolate
#ifdef DEBUG
    printf("interpolate...\n");
#endif
    int height = gm.rows, width = gm.cols;
    PERCISION grad;
    rbfres.create(height, width, CV_PERCISION);
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            grad = 0;
            for (int l = 0; l < len; ++l)
            {
                PERCISION dx = centers[l].x - w;
                PERCISION dy = centers[l].y - h;
                PERCISION dist = sqrt(dx * dx + dy * dy);
                grad += W.at<PERCISION>(l, 0) * phi(dist);
            }

            rbfres.at<PERCISION>(h, w) = grad;
        }
}

void RBF::rbf_interpolate(Mat &rbfx, Mat &rbfy, const Mat &src)
{
    vector<Point> centers;
    rbf_interpolate(rbfx, rbfy, centers, src);
}

void RBF::rbf_interpolate(Mat &rbfx, Mat &rbfy, vector<Point> &centers, const Mat &im)
{
    assert(im.channels() == 3);
    
    Mat gray;
    cvtColor(im, gray, CV_BGR2GRAY);

#ifdef DEBUG
    printf("gradient...\n");
#endif
    // gradient
    Mat gx, gy, abs_gx, abs_gy, gm;
    Sobel(gray, gx, CV_PERCISION, 1, 0);
    Sobel(gray, gy, CV_PERCISION, 0, 1);
    convertScaleAbs(gx, abs_gx);
    convertScaleAbs(gy, abs_gy);
    addWeighted(abs_gx, 0.5, abs_gy, 0.5, 0, gm, CV_PERCISION);
    
#ifdef DEBUG
    printf("calculating rbf sources...\n");
#endif
    // get source
    rbf_center(centers, gm);

#ifdef DEBUG
    printf("centers: %d\n", centers.size());
    printf("estimate weights and interpolating...\n");
#endif
    // estimate weights && interpolate
    rbf_solver(rbfx, centers, gx);
    rbf_solver(rbfy, centers, gy);
}

void RBF::plot(Mat &show, const vector<Point> &centers,
               const Mat &gx, const Mat &gy, const Mat &im, double factor)
{
    assert(im.channels() == 3);
    
    int width = im.cols, height = im.rows;
    PERCISION resolution = factor;
    PERCISION len = resolution / 2;
    int nwidth = width * resolution + len, nheight = height * resolution + len;
    show = Mat::zeros(nheight, nwidth, CV_8UC3);
    
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            PERCISION dx = gx.at<PERCISION>(h, w);
            PERCISION dy = gy.at<PERCISION>(h, w);
            PERCISION theta = atan2(dy, dx) + CV_PI / 2;
            int hh = h * resolution + len;
            int ww = w * resolution + len;
            int dw = cos(theta) * len * 2;
            int dh = sin(theta) * len * 2;
            Vec3b pixel = im.at<Vec3b>(h, w);
            Scalar color(pixel[0], pixel[1], pixel[2]);
            
            line(show, Point(ww, hh), Point(ww+dw, hh+dh), color);
            line(show, Point(ww, hh), Point(ww-dw, hh-dh), color);
        }

    for (size_t i = 0; i < centers.size(); ++i)
    {
        int x = centers[i].x * resolution + len;
        int y = centers[i].y * resolution + len;
        
        circle(show, Point(x, y), 5, Scalar(255, 0, 0), -1);
    }
}
