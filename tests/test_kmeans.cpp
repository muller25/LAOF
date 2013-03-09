#include "ML.h"
#include "Image.h"
#include "ImageIO.h"

#include <ctime>

int main(int argc, char *argv[])
{
    int height = 10, width = 2, clusters = 3, offset;
    DImage samples(width, height);
    DImage centers(width, clusters);
    UCImage labels;
    
    srand(time(NULL));
    for (int i = 0; i < height; ++i)
    {
        offset = i * width;
        for (int w = 0; w < width; ++w)
            samples[offset+w] = (int)rand() % 10;
    }

    for (int i = 0; i < clusters; ++i)
    {
        offset = i * width;
        for (int w = 0; w < width; ++w)
            centers[offset+w] = (int)rand() % 10;
    }

    printf("samples:\n");
    imprint(samples);
    
    printf("centers:\n");
    imprint(centers);

    double compactness = kmeans(centers, labels, samples, clusters);

    printf("********\n");
    printf("compactness: %.6f\n", compactness);

    printf("final centers:\n");
    imprint(centers);

    printf("labels:\n");
    imprint(labels);
    
    return 0;
}
