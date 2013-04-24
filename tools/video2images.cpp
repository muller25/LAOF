#include <highgui.h>
#include <cv.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        printf("Usage: input_video out_dir\n");
        return -1;
    }

    IplImage *frame = NULL;
    CvCapture *video = cvCreateFileCapture(argv[1]);
    
    if (video == NULL)
    {
        printf("fail to load %s\n", argv[1]);
        return -1;
    }
/*
    double fps = cvGetCaptureProperty(video, CV_CAP_PROP_FPS);
    CvSize size = cvSize(
        (int)cvGetCaptureProperty(video, CV_CAP_PROP_FRAME_WIDTH),
        (int)cvGetCaptureProperty(video, CV_CAP_PROP_FRAME_HEIGHT)
        );
    CvVideoWriter *writer = cvCreateVideoWriter(argv[2],
                                                CV_FOURCC('M', 'J', 'P', 'G'),
                                                fps,
                                                size
        );

    assert(writer != NULL);

    cvNamedWindow("video", CV_WINDOW_AUTOSIZE);
*/
    int i = 0;
    char buf[256], out[256];
    memset(out, 0, sizeof(out));
    strcat(out, argv[2]);
    strcat(out, "%03d.jpg");
    
    // load video by frames
    while (1)
    {
        frame = cvQueryFrame(video);
        if (frame == NULL)
        {
            printf("finish!\n");
            break;
        }

        sprintf(buf, out, i++);
        printf("saving in %s\n", buf);
        
        int res = cvSaveImage(buf, frame);
        if (res == 0)
        {
            printf("fail to save image\n");
            break;
        }
        
//        cvShowImage("video", frame);
//        cvWriteFrame(writer, frame);
        
        char c = cvWaitKey(33);
        if (c == 27) break;
    }

//    cvReleaseVideoWriter(&writer);
    cvReleaseCapture(&video);
//    cvDestroyWindow("video");
    
    return 0;
}
