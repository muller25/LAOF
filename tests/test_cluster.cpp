#include "ML.h"
#include "Image.h"
#include "ImageIO.h"
#include "Maths.h"

#include <ctime>

#include <cv.h>
#include <highgui.h>

int imWidth = 600, imHeight = 600;
int height = 300, width = 2;
int numOfClusters = 10;

int test_kmeans(DImage &centers, UCImage &labels, DImage &samples);
int test_SpectralCluster(DImage &centers, UCImage &labels, DImage &samples);

template <class T, class L>
void draw(Image<T> &samples, Image<L> &labels);

int main(int argc, char *argv[])
{
    DImage samples(width, height);
    DImage centers;
    UCImage labels;
    int res;
    
    randFill(samples, 10., (double)imWidth-10);
    
    // printf("samples:\n");
    // imprint(samples);

    // res = test_kmeans(centers, labels, samples);
    res = test_SpectralCluster(centers, labels, samples);

    printf("****** result ******\n");
    printf("clustesr: %d\n", res);
    
    // printf("centers:\n");
    // imprint(centers);

    // printf("labels:\n");
    // imprint(labels);

    draw(samples, labels);
    
    return 0;
}

int test_kmeans(DImage &centers, UCImage &labels, DImage &samples)
{
    kmeans(centers, labels, samples, numOfClusters);

    return numOfClusters;
}

int test_SpectralCluster(DImage &centers, UCImage &labels, DImage &samples)
{
    int n = height;
    DImage graph(n, n);
    double cost;
    
    // construct graph, or say similarity matrix
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            cost = dist2(samples.ptr()+i*width, samples.ptr()+j*width, 0, width);
            if (cost > imWidth / numOfClusters) continue;

            graph[i*n+j] = cost;
            graph[j*n+i] = graph[i*n+j];
        }
    }

    SpectralCluster(labels, graph, numOfClusters);

    centers.release();
    return numOfClusters;
}

template <class T, class L>
void draw(Image<T> &samples, Image<L> &labels)
{
    assert(samples.nWidth() == 2);
    
    const int maxNumOfClusters = 10;
    int colorIdx;
    cv::Mat img(imHeight, imWidth, CV_8UC3, cv::Scalar(255 ,255, 255));

    srand(time(NULL));
    int *color = new int[maxNumOfClusters*3];
    for (int i = 0; i < maxNumOfClusters; ++i)
    {
        color[i*3] = rand() % 25 * 10;
        color[i*3+1] = rand() % 25 * 10;
        color[i*3+2] = rand() % 25 * 10;
    }

    for (int h = 0; h < samples.nHeight(); ++h)
    {
        colorIdx = labels[h] * 3;
        
        cv::Scalar c(color[colorIdx], color[colorIdx+1], color[colorIdx+2]);
        cv::circle(img, cv::Point(samples[h*2], samples[h*2+1]), 3, c, -1);
    }

    cv::imshow("cluster result", img);
    cv::waitKey(0);
    
    delete []color;
}
