#include <stdio.h>
#include <cv.h>
#include <highgui.h>

int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        printf("Usage: image_pattern start end output_video fps\n");
        return -1;
    }

    int fps = atoi(argv[5]);
    int start = atoi(argv[2]);
    int end = atoi(argv[3]);
    char name[256];
    sprintf(name, argv[1], start);
    
    IplImage *frame = cvLoadImage(name);
    if (frame == NULL)
    {
        printf("fail to load image %s\n", name);
        return -1;
    }
    
    CvSize videoSize = cvSize(frame->width, frame->height);
    CvVideoWriter *writer = cvCreateVideoWriter(argv[4],
                                                CV_FOURCC('D', 'I', 'V', 'X'),
                                                fps,
                                                videoSize
        );

    cvReleaseImage(&frame);
 
    for (int i = start; i <= end; ++i)
    {
        sprintf(name, argv[1], i);
        printf("processing %s...\n", name);
        
        frame = cvLoadImage(name);
        if (frame == NULL)
        {
            printf("fail to load image %s\n", name);
            continue;
        }

        cvWriteFrame(writer, frame);
        cvReleaseImage(&frame);
    }

    cvReleaseVideoWriter(&writer);

    printf("done\n");
    return 0;
}
