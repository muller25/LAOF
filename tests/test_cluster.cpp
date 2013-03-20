#include "ML.h"
#include "Image.h"
#include "ImageIO.h"
#include "Maths.h"

#include <cv.h>
#include <highgui.h>

int height = 100, width = 2;

int test_kmeans2(DImage &centers, UCImage &labels, DImage &samples);
int test_SpectralCluster(DImage &centers, UCImage &labels, DImage &samples);

template <class T, class L>
void draw(Image<T> &samples, Image<L> &labels);

int main(int argc, char *argv[])
{
    DImage samples(width, height);
    DImage centers;
    UCImage labels;
    int numOfClusters;
    
    randFill(samples, 10., 290.);
    
    // printf("samples:\n");
    // imprint(samples);

    numOfClusters = test_kmeans2(centers, labels, samples);
    // numOfClusters = test_SpectralCluster(centers, labels, samples);

    printf("****** result ******\n");
    printf("clustesr: %d\n", numOfClusters);
    
    // printf("centers:\n");
    // imprint(centers);

    // printf("labels:\n");
    // imprint(labels);

    draw(samples, labels);
    
    return 0;
}

int test_kmeans2(DImage &centers, UCImage &labels, DImage &samples)
{
    return kmeans2(centers, labels, samples, 2, 10);
}

int test_SpectralCluster(DImage &centers, UCImage &labels, DImage &samples)
{
    const int numOfClusters = 3;
    int n = height;
    DImage graph(n, n);

    // construct graph, or say similarity matrix
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            graph[i*n+j] = dist2(samples.ptr()+i*width, samples.ptr()+j*width, 0, width);
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
    cv::Mat img(300, 300, CV_8UC3, cv::Scalar(255 ,255, 255));

    int *color = new int[maxNumOfClusters];
    int interval = 256 / maxNumOfClusters;
    for (int i = 0; i < maxNumOfClusters; ++i)
        color[i] = i * interval;

    for (int h = 0; h < samples.nHeight(); ++h)
    {
        colorIdx = color[(int)labels[h*2]];
        
        cv::Scalar c(colorIdx, colorIdx, colorIdx);
        cv::circle(img, cv::Point(samples[h*2], samples[h*2+1]), 2, c, -1);
    }

    cv::imshow("cluster result", img);
    cv::waitKey(0);
    
    delete []color;
}
