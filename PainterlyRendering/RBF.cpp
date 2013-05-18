#include "RBF.h"

#include <stdio.h>

#include <highgui.h>
using namespace cv;

const double PI = atan(1.) * 4;

double RBF::factor = 0.125;
int RBF::wsize = std::max(5 * RBF::factor, 1.);
int RBF::range = std::max(10 * RBF::factor, 1.);

/*
  描述：
  计算RBF的中心点，中心点必须符合
  1. 梯度幅度大于等于阈值
  2. 必须为局部最大值
  3. 相邻中心点的距离大于等于最小距离
  参数：
  centers: 返回值，每个元素由 (x, y) 组成
  gm: 梯度幅度
  thres: 阈值，当梯度幅度大于等于阈值时才认为该点为候选中心点
  wsize: 极大值窗口
  range： 中心点之间的最小距离
*/
void RBF::rbf_center(vector<Point> &centers, const Mat &gm,
                     double thres, int wsize, int range)
{
    centers.clear();

    double minVal, maxVal;
    minMaxLoc(gm, &minVal, &maxVal);
    printf("magnitude: %.6f .. %.6f\n", minVal, maxVal);
    
    printf("finding centers...\n");
    int width = gm.cols, height = gm.rows;
    bool maxima, crowd;
    Mat taken = Mat::zeros(height, width, CV_8U);
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            double mag = fabs(gm.at<PERCISION>(h, w));
            if (mag < thres) continue;

            // 极大值测试
            maxima = true;
            for (int hh = -wsize; hh <= wsize && maxima; ++hh)
                for (int ww = -wsize; ww <= wsize; ++ww)
                {
                    int nh = h + hh;
                    int nw = w + ww;

                    nh = std::min(std::max(nh, 0), height-1);
                    nw = std::min(std::max(nw, 0), width-1);
                    if (mag < gm.at<PERCISION>(nh, nw)){
                        maxima = false;
                        break;
                    }
                }

            if (!maxima) continue;

            // 中心点距离测试
            crowd = false;
            for (int hh = -range; hh <= range && !crowd; ++hh)
                for (int ww = -range; ww <= range; ++ww)
                {
                    int nh = h + hh;
                    int nw = w + ww;

                    nh = std::min(std::max(nh, 0), height-1);
                    nw = std::min(std::max(nw, 0), width-1);
                    if (taken.at<uchar>(nh, nw) == 1)
                    {
                        crowd = true;
                        break;
                    }
                }

            if (crowd) continue;

            centers.push_back(Point(w, h));
            taken.at<uchar>(h, w) = 1;
        }
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
    printf("interpolate...\n");
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

void RBF::rbf(Mat &res, const Mat &src, int radius)
{
    assert(src.channels() == 3);
    
    Mat im, gray;
    resize(src, im, Size(0, 0), factor, factor);
    cvtColor(im, gray, CV_BGR2GRAY);
    
    // smooth
    printf("smoothing...\n");
    Mat smooth;
    int ksize = 2 * radius + 1;
    GaussianBlur(gray, smooth, Size(ksize, ksize), radius, radius);

    // gradient
    printf("gradient...\n");
    Mat gx, gy, gm;
    Sobel(smooth, gx, CV_PERCISION, 1, 0);
    Sobel(smooth, gy, CV_PERCISION, 0, 1);
    magnitude(gx, gy, gm);
    
    // get source
    printf("calculating rbf sources...\n");
    int width = im.cols, height = im.rows, len;
    vector<Point> centersx, centersy, centers;

    rbf_center(centers, gm, 100, wsize, range);
    len = centers.size();
    printf("centers: %d\n", len);
    
    // estimate weights && interpolate
    printf("estimate weights and interpolating...\n");
    Mat rbfx, rbfy;
    rbf_solver(rbfx, centers, gx);
    rbf_solver(rbfy, centers, gy);
    res = Mat::zeros(height, width, CV_MAKETYPE(CV_PERCISION, 2));
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            Vec<PERCISION, 2> grad;
            grad[0] = rbfx.at<PERCISION>(h, w);
            grad[1] = rbfy.at<PERCISION>(h, w);
            res.at<Vec<PERCISION, 2> >(h, w) = grad;
        }
    }

    Mat origin, rbfres;
    plot_gradient(origin, gx, gy, im, factor);
    plot_gradient(rbfres, rbfx, rbfy, im, factor);
    imshow("origin", origin);
    imshow("rbf", rbfres);
    waitKey(0);
}

void RBF::orientation(Mat &field, const Mat &src, int radius)
{
    Mat rbfres;
    rbf(rbfres, src, radius);
    orientation(field, rbfres);
}

void RBF::orientation(Mat &field, const Mat &rbfres)
{
    assert(rbfres.channels() == 2);
    
    int width = rbfres.cols, height = rbfres.rows;
    field = Mat::zeros(height, width, CV_PERCISION);
    
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            Vec<PERCISION, 2> gradient = rbfres.at<Vec<PERCISION, 2> >(h, w);
            field.at<PERCISION>(h, w) = atan2(gradient[1], gradient[0]) + PI / 2;
        }
}

void RBF::plot_gradient(Mat &show, const Mat &gx, const Mat &gy,
                   const Mat &im, double factor)
{
    assert(im.channels() == 3);
    
    int width = im.cols, height = im.rows;
    PERCISION resolution = 1 / factor;
    PERCISION len = resolution / 2;
    int nwidth = width * resolution + len, nheight = height * resolution + len;
    show = Mat::zeros(nheight, nwidth, CV_8UC3);
    
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            PERCISION dx = gx.at<PERCISION>(h, w);
            PERCISION dy = gy.at<PERCISION>(h, w);
            PERCISION theta = atan2(dy, dx) + PI / 2;
            int hh = h * resolution + len;
            int ww = w * resolution + len;
            int dw = cos(theta) * len * 2;
            int dh = sin(theta) * len * 2;
            Vec3b pixel = im.at<Vec3b>(h, w);
            Scalar color(pixel[0], pixel[1], pixel[2]);
            
            line(show, Point(ww, hh), Point(ww+dw, hh+dh), color);
            line(show, Point(ww, hh), Point(ww-dw, hh-dh), color);
        }
}
