#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <cstdio>
#include <vector>
using std::vector;

typedef Vec2i POSITION;
typedef Vec2f MOTION;
typedef Vec3b PIXEL;

// #define DEBUG

bool imreadf(Mat &m, const char *filename)
{
    cv::FileStorage fs(filename, cv::FileStorage::READ);

    if (!fs.isOpened())
    {
        fs.release();
        return false;
    }

    fs["matrix"] >> m;
    fs.release();
    
    return true;
}

void covariance(Mat &cov, const Mat &x, const Mat &y)
{
    Mat arr[] = {x, y};
    Mat mean;
    calcCovarMatrix(arr, 2, cov, mean, CV_COVAR_NORMAL, CV_64F);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("./trajectory inIm inFile outDir\n");
        return -1;
    }
    
    char inIm[256], inFile[256], outDir[256] = {0}, buf[256];
    strcpy(inIm, argv[1]);
    strcpy(inFile, argv[2]);
    strcpy(outDir, argv[3]);
    strcat(outDir, "%s%03d.jpg");
    
    const int timeWin = 24;
    const int maxDiff = 50 * 3;
    const double threshold = 1e-3;
    
    Mat u, v;
    sprintf(buf, inFile, "u", 0);
    imreadf(u, buf);
    sprintf(buf, inFile, "v", 0);
    imreadf(v, buf);

    int width = u.cols, height = u.rows;
    Mat trust = Mat::ones(height, width, CV_8U);
    Mat curPos(height, width, CV_32SC2); // current position of trajectory
    Mat trajectory(height * width, timeWin, CV_32FC2);

    // init trajectory
    printf("init trajectory...\n");
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            curPos.at<POSITION>(h, w) = POSITION(w, h);

            double uVal = cvRound(u.at<double>(h, w));
            double vVal = cvRound(v.at<double>(h, w));
            trajectory.at<MOTION>(h*width+w, 0) = MOTION(uVal, vVal);
        }

    // calculate trajectory
    Mat curIm, nextIm, cov, failWarp;

    printf("calculate trajectory...\n");
    for (int t = 1; t < timeWin; ++t)
    {
        sprintf(buf, inIm, t-1);
        curIm = imread(buf);
        sprintf(buf, inIm, t);
        nextIm = imread(buf);
        sprintf(buf, inFile, "u", t);
        imreadf(u, buf);
        sprintf(buf, inFile, "v", t);
        imreadf(v, buf);

        covariance(cov, u, v);

#ifdef DEBUG
        double minVal, maxVal;
        minMaxLoc(cov, &minVal, &maxVal);
        printf("cov: %.6f .. %.6f\n", minVal, maxVal);
        sprintf(buf, "cov", t);
        imwrite(buf, (cov - minVal) / (maxVal - minVal));

        sprintf(buf, "trust", t);
        imwrite(buf, (cov <= threshold));

        failWarp = Mat::zeros(height, width, CV_8U);
#endif
        
        for (int h = 0; h < height; ++h)
            for (int w = 0; w < width; ++w)
            {
                int offset = h * width + w;
                POSITION pos = curPos.at<POSITION>(h, w);
                MOTION motion = trajectory.at<MOTION>(offset, t-1);
                PIXEL color = curIm.at<PIXEL>(h, w);
                int nw = pos[0] + motion[0];
                int nh = pos[1] + motion[1];
                
                if (nw < 0 || nw >= width || nh >= height || nh < 0)
                {
                    trust.at<uchar>(h, w) = 0;
                    trajectory.at<MOTION>(offset, t) = motion;
                    int diff = color[0] + color[1] + color[2];
                    if (diff > maxDiff) {
                        trust.at<uchar>(h, w) = 0;
                        failWarp.at<uchar>(h, w) = 1;
                    }
                }
                else {
                    if (cov.at<double>(h, w) > threshold) trust.at<uchar>(h, w) = 0;

                    PIXEL ncolor = nextIm.at<PIXEL>(nh, nw);
                    int diff = abs(ncolor[0] - color[0]) + abs(ncolor[1] - color[1]) + abs(ncolor[2] - color[2]);
                    if (diff > maxDiff) {
                        trust.at<uchar>(h, w) = 0;
                        failWarp.at<uchar>(h, w) = 1;
                    }
                    
                    curPos.at<POSITION>(h, w) = POSITION(nw, nh);
                    double uVal = cvRound(u.at<double>(nh, nw));
                    double vVal = cvRound(v.at<double>(nh, nw));
                    trajectory.at<MOTION>(offset, t) = MOTION(uVal, vVal);
                }
            }
    }

    printf("sampling confident trajectory...\n");
    int idx = 0, nsamples = countNonZero(trust);
    vector<int> id2tid(width*height, -1);
    Mat samples(nsamples, timeWin, CV_32FC2);
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            if (trust.at<uchar>(h, w) == 0) continue;

            int offset = h * width + w;
            for (int t = 0; t < timeWin; ++t)
                samples.at<MOTION>(idx, t) = trajectory.at<MOTION>(offset, t);

            id2tid[offset] = idx++;
        }
    printf("trust samples: %d / %d\n", nsamples, height*width);
    
    Mat label;
    const int nlabel = 3;
    printf("running kmeans....\n");
    kmeans(samples, nlabel, label, TermCriteria(TermCriteria::COUNT + TermCriteria::EPS, 100, 1), 5, KMEANS_PP_CENTERS);
    curIm = Mat::zeros(height, width, CV_8U);
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            int offset = h * width + w;
            int tid = id2tid[offset];
            if (tid >= 0)
                curIm.at<uchar>(h, w) = (label.at<int>(tid) + 1) * 255 / (nlabel+1);
        }
    
    sprintf(buf, outDir, "label", 0);
    imwrite(buf, curIm);
    sprintf(buf, outDir, "trust", 0);
    imwrite(buf, trust);
    
    for (int t = 0; t < timeWin; ++t)
    {
        nextIm = Mat::zeros(height, width, CV_8U);
        for (int h = 0; h < height; ++h)
            for (int w = 0; w < width; ++w)
            {
                int offset = h * width + w;
                MOTION motion = trajectory.at<MOTION>(offset, t);
                int nw = w + motion[0], nh = h + motion[1];
                if (nw < 0 || nw >= width || nh >= height || nh < 0) continue;

                nextIm.at<uchar>(nh, nw) = curIm.at<uchar>(h, w);
            }

        sprintf(buf, outDir, "label", t+1);
        imwrite(buf, nextIm);
        nextIm.copyTo(curIm);
    }
    
    return 0;
}
