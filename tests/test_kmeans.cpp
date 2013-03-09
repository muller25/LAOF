#include "ML.h"
#include "Image.h"
#include "ImageIO.h"

#include <ctime>

int main(int argc, char *argv[])
{
    int height = 20, width = 2, clusters = 10, offset;
    DImage samples(width, height);
    DImage centers;
    UCImage labels;
    
    srand(time(NULL));
    for (int i = 0; i < height; ++i)
    {
        offset = i * width;
        for (int w = 0; w < width; ++w)
            samples[offset+w] = rand() % height;
    }

    printf("samples:\n");
    imprint(samples);
    
    // printf("centers:\n");
    // imprint(centers);

    // double compactness = kmeans(centers, labels, samples, clusters);
    clusters = kmeans2(centers, labels, samples, 2, 20);
    
    printf("********\n");
    printf("final clustesr: %d\n", clusters);
    
    // printf("compactness: %.6f\n", compactness);

    printf("final centers:\n");
    imprint(centers);

    printf("labels:\n");
    imprint(labels);
    
    return 0;
}
