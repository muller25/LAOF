#include <stdio.h>
#include <cv.h>
#include <highgui.h>

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printf("Usage: sub-video-pattern start end output-video\n");
        return -1;
    }

    IplImage *frame;
    char name[256];
    int start = atoi(argv[2]);
    int end = atoi(argv[3]);
    sprintf(name, argv[1], start);
    CvCapture *video = cvCreateFileCapture(name);    
    if (video == NULL)
    {
        printf("fail to load video %s\n", name);
        return -1;
    }


    CvSize videoSize = cvSize((int)cvGetCaptureProperty(video, CV_CAP_PROP_FRAME_WIDTH),
                              (int)cvGetCaptureProperty(video, CV_CAP_PROP_FRAME_HEIGHT)
        );

    double fps = cvGetCaptureProperty(video, CV_CAP_PROP_FPS);
    cvReleaseCapture(&video);    

    CvVideoWriter *writer = cvCreateVideoWriter(argv[argc - 1],
                                                CV_FOURCC('D', 'I', 'V', 'X'),
                                                fps,
                                                videoSize
        );

    if (writer == NULL)
    {
        printf("fail to create writer for video %s\n", argv[argc - 1]);
        return -1;
    }


    for (int i = start; i <= end; ++i)
    {
        sprintf(name, argv[1], i);
        printf("processing %s...\n", name);
        video = cvCreateFileCapture(name);

        if (video == NULL)
        {
            printf("fail to load video %s\n", name);
            return -1;
        }
        
        while (frame = cvQueryFrame(video), frame != NULL)
            cvWriteFrame(writer, frame);

        cvReleaseCapture(&video);
    }

    cvReleaseVideoWriter(&writer);

    printf("done!\n");
    return 0;
}
